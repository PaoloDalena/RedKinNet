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
#' A component is matched to the first rule in \code{rules_db} whose size constraint
#' equals \code{nrow(fam_df)} and whose \code{match_fn(fam_df)} evaluates to \code{TRUE}.
#' If no rule matches, all nodes are returned as discarded.
#'
#' Randomization uses one of two mechanisms:
#' \enumerate{
#'   \item If \code{rule$alpha > 1} and \code{g_sub} is provided, the function samples
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
#' @param rules_db A list of rule objects. Each rule is expected to provide:
#'   \itemize{
#'     \item \code{name}: rule name.
#'     \item \code{size_allowed}: integer size constraint.
#'     \item \code{match_fn}: function(\code{fam_df}) returning a scalar logical.
#'     \item \code{kept_rules}: list of character expressions.
#'     \item \code{discarded_rules}: list of character expressions.
#'     \item \code{conflict_rules}: list of character vectors, each defining a candidate group.
#'     \item Optional: \code{alpha} (numeric), \code{kept_fn}, \code{discarded_fn}
#'           for graph-dependent deterministic selection.
#'   }
#' @param seed Optional integer RNG seed used for randomized choices.
#' @param return_diagnostics Logical; if \code{TRUE}, returns diagnostics in addition to
#'   \code{kept} and \code{discarded}.
#' @param g_sub Optional \pkg{igraph} graph for the same component. Required if using
#'   the maximum independent set mechanism (\code{rule$alpha > 1}). Vertex names
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
#' out <- rkn_select_random_upto5(fam_df, rules_db, seed = 1)
#' out_diag <- rkn_select_random_upto5(fam_df, rules_db, seed = 1,
#'                                    return_diagnostics = TRUE, g_sub = g_sub)
#' }
#'
#' @export
rkn_select_random_upto5 <- function(fam_df,
                                    rules_db,
                                    seed = NULL,
                                    return_diagnostics = FALSE,
                                    g_sub = NULL) {

  # Harmonize identifier columns.
  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  # Match the component to a rule.
  rule <- find_matching_rule(fam_df, rules_db)

  # Build per-node diagnostics (upto5; aligned with 6plus fields).
  build_node_diag <- function(ids,
                              kept_apriori, discarded_apriori,
                              candidate_simple, candidate_complex,
                              kept_final, discarded_final,
                              pre_reason_map = NULL) {

    pre_class <- rep(NA_character_, length(ids))
    names(pre_class) <- ids

    pre_class[ids %in% kept_apriori]      <- "kept_apriori"
    pre_class[ids %in% discarded_apriori] <- "discarded_apriori"
    pre_class[ids %in% candidate_simple]  <- "candidate_simple"
    pre_class[ids %in% candidate_complex] <- "candidate_complex"

    # Fallback (should be rare): any remaining node is treated as complex candidate.
    pre_class[is.na(pre_class)] <- "candidate_complex"

    post_class <- rep(NA_character_, length(ids))
    names(post_class) <- ids
    post_class[intersect(ids, kept_final)]      <- "kept"
    post_class[intersect(ids, discarded_final)] <- "discarded"
    post_class[is.na(post_class)] <- "discarded"

    pre_reason <- rep(NA_character_, length(ids))
    names(pre_reason) <- ids
    if (!is.null(pre_reason_map) && length(pre_reason_map) > 0) {
      pre_reason[names(pre_reason_map)] <- unname(pre_reason_map)
    }

    data.frame(
      id = ids,
      pre_class = unname(pre_class[ids]),
      post_class = unname(post_class[ids]),
      is_selected = unname(post_class[ids] == "kept"),
      pre_reason = unname(pre_reason[ids]),
      stringsAsFactors = FALSE
    )
  }

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

  # Deterministic (table-driven) selection.
  kept_apriori      <- apply_rule_filter(fam_df, unlist(rule$kept_rules))
  discarded_apriori <- apply_rule_filter(fam_df, unlist(rule$discarded_rules))

  # Optional graph-dependent deterministic selection (kept_fn / discarded_fn).
  if (!is.null(g_sub)) {
    if (!is.null(rule$kept_fn) && is.function(rule$kept_fn)) {
      kept_graph <- rule$kept_fn(fam_df, g_sub)
      kept_apriori <- unique(c(kept_apriori, kept_graph))
    }
    if (!is.null(rule$discarded_fn) && is.function(rule$discarded_fn)) {
      disc_graph <- rule$discarded_fn(fam_df, g_sub)
      discarded_apriori <- unique(c(discarded_apriori, disc_graph))
    }
  }

  # Conflict groups (always excluding a priori decisions).
  conflict_groups <- list()
  candidates_all  <- character(0)

  if (length(rule$conflict_rules) > 0) {
    for (grp in rule$conflict_rules) {
      cand_nodes <- apply_rule_filter(fam_df, grp)
      cand_nodes <- setdiff(cand_nodes, union(kept_apriori, discarded_apriori))
      cand_nodes <- unique(cand_nodes)
      if (length(cand_nodes) > 0) conflict_groups[[length(conflict_groups) + 1]] <- cand_nodes
    }
    if (length(conflict_groups) > 0) {
      candidates_all <- unique(unlist(conflict_groups))
    }
  }

  # Randomized resolution.
  used_largest_ivs <- FALSE
  kept_afterrandomization <- character(0)
  discarded_afterrandomization <- character(0)

  if (!is.null(rule$alpha) && isTRUE(rule$alpha > 1) && !is.null(g_sub)) {

    # Eligible nodes: all nodes except discarded a priori; optionally exclude kept a priori
    # and their 1-hop neighbors.
    eligible <- setdiff(fam_df$.name, discarded_apriori)

    if (length(kept_apriori) > 0) {
      nb_idx  <- igraph::neighborhood(g_sub, order = 1, nodes = kept_apriori)
      nb_names <- unique(setdiff(
        unlist(lapply(nb_idx, function(vs) igraph::V(g_sub)$name[vs])),
        kept_apriori
      ))
      eligible <- setdiff(eligible, union(kept_apriori, nb_names))
    }

    g_ind <- igraph::induced_subgraph(g_sub, vids = eligible)
    sets  <- igraph::largest_ivs(g_ind)
    set_names <- lapply(sets, function(s) igraph::V(g_ind)$name[s])

    if (!is.null(seed)) set.seed(seed)
    chosen_names <- if (length(set_names) > 0) set_names[[sample(seq_along(set_names), 1)]] else character(0)

    kept_afterrandomization <- unique(chosen_names)

    # If conflict groups were not defined, use eligible nodes for diagnostics.
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

  # Final reconciliation (kept overrides discarded).
  kept_final <- unique(c(kept_apriori, kept_afterrandomization))
  discarded_final <- unique(c(discarded_apriori, discarded_afterrandomization))
  discarded_final <- setdiff(discarded_final, kept_final)

  # Candidate complexity flag (complex if multiple groups or alpha > 1).
  is_complex <- (length(conflict_groups) > 1) ||
    (!is.null(rule$alpha) && isTRUE(rule$alpha > 1))

  candidate_complex <- if (is_complex) unique(candidates_all) else character(0)
  candidate_simple  <- if (!is_complex) unique(candidates_all) else character(0)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  }

  # Minimal pre-reason map for diagnostics/plotting.
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
#' Expressions are evaluated via \code{.node_filter()}, which treats columns of
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
apply_rule_filter <- function(fam_df, rule_exprs) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (is.null(fam_df$.name)) stop("`fam_df` must contain a `.name` column.")
  if (length(rule_exprs) == 0) return(character(0))

  out <- character(0)

  for (expr_txt in rule_exprs) {
    mask <- .node_filter(fam_df, expr_txt)

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
find_matching_rule <- function(fam_df, rules_db) {
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
.node_filter <- function(fam_df, expr_txt) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (!is.character(expr_txt) || length(expr_txt) != 1L) {
    stop("`expr_txt` must be a single character string.")
  }

  if (identical(expr_txt, "TRUE"))  return(rep(TRUE,  nrow(fam_df)))
  if (identical(expr_txt, "FALSE")) return(rep(FALSE, nrow(fam_df)))

  eval(parse(text = expr_txt), envir = as.list(fam_df))
}
