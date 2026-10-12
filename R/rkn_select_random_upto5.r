#' Random selection within connected components of size <= 5
#'
#' @description
#' Applies the RedKinNet randomized node selection strategy to a connected
#' component with at most 5 nodes. The function combines:
#' \itemize{
#'   \item deterministic (a priori) keep/discard decisions based on a matched rule; and
#'   \item randomized resolution of remaining ambiguous conflict groups.
#' }
#'
#' A component is matched to the first rule in \code{rules_db} (or in \code{custom_rules_db} if provided) whose size constraint
#' equals \code{nrow(fam_df)} and whose \code{match_fn(fam_df)} evaluates to \code{TRUE}.
#' If no rule matches, all nodes are returned as discarded.
#'
#' Randomization uses one of two mechanisms:
#' \enumerate{
#'   \item If \code{rule$use_largest_ivs} is \code{TRUE} and \code{g_sub} is provided, the function samples
#'   uniformly among maximum independent vertex sets returned by
#'   \code{igraph::largest_ivs()} on an induced subgraph of eligible vertices.
#'   \item Otherwise, it performs simple random sampling within each conflict group
#'   defined by \code{rule$conflict_rules}.
#' }
#'
#' When \code{return_diagnostics = TRUE}, the function returns per-node and per-family
#' diagnostics aligned with the 6+ selector output (fields \code{candidate_simple},
#' \code{candidate_complex}, \code{node_diagnostics}, \code{family_diagnostics}).
#'
#' @param fam_df A data frame describing nodes in the component. Must contain a node
#'   identifier column \code{.name} (preferred) or \code{name} (promoted to \code{.name}).
#'   Must also contain any columns referenced by rule expressions in \code{rules_db}.
#' @param custom_rules_db Optional list of rule objects. If provided, it overrides the
#'   internal rule database used for components of size <= 5. If \code{NULL} (default),
#'   the function uses the package dataset \code{rules_db}; see \code{\link{rules_db}}.
#'
#'   Each rule is expected to provide:
#'   \itemize{
#'     \item \code{name}: rule name.
#'     \item \code{size_allowed}: integer size constraint.
#'     \item \code{match_fn}: function(\code{fam_df}) returning a scalar logical.
#'     \item \code{kept_rules}: list of character expressions.
#'     \item \code{discarded_rules}: list of character expressions.
#'     \item \code{conflict_rules}: list of character vectors, each defining a candidate group.
#'     \item Optional: \code{use_largest_ivs} (logical), \code{kept_fn}, \code{discarded_fn}
#'           for graph-dependent deterministic selection.
#'   }
#' @param seed Optional integer RNG seed used for randomized choices.
#' @param return_diagnostics Logical; if \code{TRUE}, returns diagnostics in addition to
#'   \code{kept} and \code{discarded}.
#' @param g_sub Optional \pkg{igraph} graph for the same component. Required if using
#'   the maximum independent set mechanism (\code{rule$use_largest_ivs}). Vertex names
#'   \code{igraph::V(g_sub)$name} must match \code{fam_df$.name}.
#'
#' @return
#' If \code{return_diagnostics = FALSE}, a list with:
#' \itemize{
#'   \item \code{kept}: character vector of selected node names.
#'   \item \code{discarded}: character vector of discarded node names.
#' }
#'
#' If \code{return_diagnostics = TRUE}, the list also includes:
#' \itemize{
#'   \item \code{structure_name}: matched rule name (or \code{"unknown"}).
#'   \item \code{kept_apriori}, \code{discarded_apriori}.
#'   \item \code{candidate_simple}, \code{candidate_complex}.
#'   \item \code{kept_afterrandomization}, \code{discarded_afterrandomization}.
#'   \item \code{used_largest_ivs}.
#'   \item \code{node_diagnostics}: one row per node (pre/post class and reason).
#'   \item \code{family_diagnostics}: one-row summary data frame.
#' }
#'
#' @examples
#' \dontrun{
#' out <- rkn_select_random_upto5(fam_df, seed = 1)
#' out_diag <- rkn_select_random_upto5(fam_df, seed = 1,
#'                                    return_diagnostics = TRUE, g_sub = g_sub)
#' }
#'
#' @export
rkn_select_random_upto5 <- function(fam_df,
                                    custom_rules_db = NULL,
                                    seed = NULL,
                                    return_diagnostics = FALSE,
                                    g_sub = NULL) {

  if (!is.null(seed)) withr::local_preserve_seed()

  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  rules_db_use <- if (is.null(custom_rules_db)) RedKinNet::rules_db else custom_rules_db
  rule <- .rkn_find_matching_rule(fam_df, rules_db_use)
  build_node_diag <- .rkn_build_node_diag_upto5

  # No matching rule: discard all nodes.
  if (is.null(rule)) {
    kept_final <- character(0)
    discarded_final <- fam_df$.name

    if (!return_diagnostics) {
      return(list(kept = kept_final, discarded = discarded_final))
    }

    node_diagnostics <- data.frame(
      id = fam_df$.name,
      pre_class = rep("discarded_apriori", nrow(fam_df)),
      post_class = rep("discarded", nrow(fam_df)),
      is_selected = rep(FALSE, nrow(fam_df)),
      pre_reason = rep("no_matching_rule", nrow(fam_df)),
      stringsAsFactors = FALSE
    )

    family_diagnostics <- data.frame(
      n_nodes = nrow(fam_df),
      used_largest_ivs = FALSE,
      structure_name = "unknown",
      stringsAsFactors = FALSE
    )

    return(list(
      structure_name = "unknown",
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = FALSE,
      node_diagnostics = node_diagnostics,
      family_diagnostics = family_diagnostics
    ))
  }

  ac <- .rkn_apriori_and_conflicts_upto5(fam_df, rule, g_sub)
  kept_apriori      <- ac$kept_apriori
  discarded_apriori <- ac$discarded_apriori
  conflict_groups   <- ac$conflict_groups
  candidates_all    <- ac$candidates_all

  used_largest_ivs <- FALSE
  kept_afterrandomization <- character(0)
  discarded_afterrandomization <- character(0)

  if (isTRUE(rule$use_largest_ivs) && !is.null(g_sub)) {

    # Nodes adjacent to an a-priori kept node can no longer be selected.
    eligible <- setdiff(fam_df$.name, discarded_apriori)
    eligible <- .rkn_exclude_neighbors_of_kept(g_sub, eligible, kept_apriori)

    g_ind <- igraph::induced_subgraph(g_sub, vids = eligible)
    sets  <- igraph::largest_ivs(g_ind)
    ind_names <- igraph::V(g_ind)$name
    set_names <- vector("list", length(sets))
    for (j in seq_along(sets)) set_names[[j]] <- ind_names[sets[[j]]]

    if (!is.null(seed)) set.seed(seed)
    chosen_names <- if (length(set_names) > 0) set_names[[sample(seq_along(set_names), 1)]] else character(0)

    kept_afterrandomization <- unique(chosen_names)

    if (length(candidates_all) == 0) candidates_all <- eligible

    discarded_afterrandomization <- setdiff(candidates_all, kept_afterrandomization)
    used_largest_ivs <- TRUE

  } else if (length(conflict_groups) > 0) {

    if (!is.null(seed)) set.seed(seed)

    if (length(conflict_groups) == 1) {
      grp <- conflict_groups[[1]]
      kept_afterrandomization <- if (length(grp) == 1) grp else sample(grp, size = 1)
    } else {
      U <- unique(unlist(conflict_groups))
      k <- length(conflict_groups)

      if (length(U) >= k) {
        kept_afterrandomization <- sample(U, size = k, replace = FALSE)
      } else {
        chosen <- character(0)
        for (grp in conflict_groups) {
          avail <- setdiff(grp, chosen)
          if (length(avail) > 0) chosen <- c(chosen, sample(avail, 1))
        }
        kept_afterrandomization <- unique(chosen)
      }
    }

    discarded_afterrandomization <- setdiff(candidates_all, kept_afterrandomization)
  }

  kept_afterrandomization <- unique(kept_afterrandomization)

  # Kept overrides discarded.
  kept_final <- unique(c(kept_apriori, kept_afterrandomization))
  discarded_final <- unique(c(discarded_apriori, discarded_afterrandomization))
  discarded_final <- setdiff(discarded_final, kept_final)

  is_complex <- (length(conflict_groups) > 1) || isTRUE(rule$use_largest_ivs)

  candidate_complex <- if (is_complex) unique(candidates_all) else character(0)
  candidate_simple  <- if (!is_complex) unique(candidates_all) else character(0)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  }

  pre_reason_map <- character(0)
  if (used_largest_ivs && length(candidate_complex) > 0) {
    pre_reason_map[candidate_complex] <- "largest_ivs_fallback"
  } else if (length(conflict_groups) > 1 && length(candidate_complex) > 0) {
    pre_reason_map[candidate_complex] <- "multi_conflict_groups"
  } else if (length(conflict_groups) == 1 && length(candidate_simple) > 0) {
    pre_reason_map[candidate_simple] <- "single_conflict_group"
  }

  node_diagnostics <- build_node_diag(
    ids = fam_df$.name,
    kept_apriori = unique(kept_apriori),
    discarded_apriori = unique(discarded_apriori),
    candidate_simple = unique(candidate_simple),
    candidate_complex = unique(candidate_complex),
    kept_final = kept_final,
    discarded_final = discarded_final,
    pre_reason_map = pre_reason_map
  )

  family_diagnostics <- data.frame(
    n_nodes = nrow(fam_df),
    used_largest_ivs = used_largest_ivs,
    structure_name = rule$name,
    stringsAsFactors = FALSE
  )

  list(
    structure_name = rule$name,
    kept_apriori = unique(kept_apriori),
    discarded_apriori = unique(discarded_apriori),
    candidate_simple = unique(candidate_simple),
    candidate_complex = unique(candidate_complex),
    kept_afterrandomization = unique(kept_afterrandomization),
    discarded_afterrandomization = unique(discarded_afterrandomization),
    kept = kept_final,
    discarded = discarded_final,
    used_largest_ivs = used_largest_ivs,
    node_diagnostics = node_diagnostics,
    family_diagnostics = family_diagnostics
  )
}

#' Internal helper: apply a set of rule expressions to a component data frame
#'
#' @description
#' Evaluates one or more logical expressions (as text) in the context of \code{fam_df}
#' and returns the corresponding node names (from \code{fam_df$.name}) for which the
#' expression evaluates to \code{TRUE}.
#'
#' Expressions are evaluated via \code{.rkn_node_filter()}, which treats columns of
#' \code{fam_df} as vectors in the evaluation environment.
#'
#' @param fam_df A data frame describing a connected component. Must contain a
#'   character column \code{.name}.
#' @param rule_exprs Character vector of expressions. Each element must evaluate
#'   to a logical vector of length \code{nrow(fam_df)}, or to a scalar logical
#'   (recycled to length \code{nrow(fam_df)}).
#'
#' @return A character vector of unique node names matching at least one expression.
#'
#' @keywords internal
#' @noRd
.rkn_apply_rule_filter <- function(fam_df, rule_exprs) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (is.null(fam_df$.name)) stop("`fam_df` must contain a `.name` column.")
  if (length(rule_exprs) == 0) return(character(0))

  out <- character(0)

  for (expr_txt in rule_exprs) {
    mask <- .rkn_node_filter(fam_df, expr_txt)

    if (length(mask) == 1L) {
      mask <- rep(isTRUE(mask), nrow(fam_df))
    } else {
      mask <- as.logical(mask)
      mask[is.na(mask)] <- FALSE
    }

    out <- c(out, fam_df$.name[which(mask)])
  }

  unique(out)
}

#' Internal helper: find the first rule matching a component
#'
#' @description
#' Scans \code{rules_db} and returns the first rule whose \code{size_allowed}
#' equals \code{nrow(fam_df)} and whose \code{match_fn(fam_df)} returns \code{TRUE}.
#'
#' @param fam_df A data frame describing a connected component.
#' @param rules_db A list of rule objects.
#'
#' @return The first matching rule (a list), or \code{NULL} if no rule matches.
#'
#' @keywords internal
#' @noRd
.rkn_find_matching_rule <- function(fam_df, rules_db) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (!is.list(rules_db)) stop("`rules_db` must be a list.")

  for (rule in rules_db) {
    if (isTRUE(rule$size_allowed == nrow(fam_df)) && isTRUE(rule$match_fn(fam_df))) {
      return(rule)
    }
  }

  NULL
}

#' Internal helper: evaluate a rule expression in the context of a component
#'
#' @description
#' Evaluates \code{expr_txt} in an environment where columns of \code{fam_df}
#' are exposed as vectors.
#'
#' @param fam_df A data frame describing a connected component.
#' @param expr_txt A single character string containing an R expression.
#'
#' @return The raw evaluation result (typically a logical vector). The caller is
#'   responsible for coercion and NA handling.
#'
#' @details
#' This helper uses \code{eval(parse())}. It is intended for controlled rule inputs
#' coming from an internal rule database, not for untrusted user strings.
#'
#' @keywords internal
#' @noRd
.rkn_node_filter <- function(fam_df, expr_txt) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (!is.character(expr_txt) || length(expr_txt) != 1L) {
    stop("`expr_txt` must be a single character string.")
  }

  if (identical(expr_txt, "TRUE"))  return(rep(TRUE,  nrow(fam_df)))
  if (identical(expr_txt, "FALSE")) return(rep(FALSE, nrow(fam_df)))

  eval(parse(text = expr_txt), envir = as.list(fam_df))
}
