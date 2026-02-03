#' Random selection within connected components of size <= 5
#'
#' @description
#' Selects a subset of nodes to keep (and optionally discard) within a connected
#' component of size \eqn{\le 5}, combining:
#' (i) deterministic (a priori) rules, and
#' (ii) randomized resolution of remaining conflicts.
#'
#' The function matches the component against a rule in \code{rules_db}. If a rule
#' is found, it applies deterministic keep/discard filters and then resolves any
#' residual ambiguity either by:
#' \itemize{
#'   \item choosing a random maximum independent vertex set when \code{rule$alpha > 1}
#'         and \code{g_sub} is provided, or
#'   \item sampling one node per conflict group (or an equivalent randomized fallback).
#' }
#'
#' If no rule matches, all nodes are returned as discarded (and none kept).
#'
#' @param fam_df A data frame describing the nodes in the connected component.
#'   Must include a node identifier column \code{.name} (preferred) or \code{name}
#'   (the function will align them), and any columns referenced by \code{rules_db}
#'   expressions (e.g., \code{degree}, \code{max_path_length}).
#'
#' @param rules_db A list of rule objects. Each rule is expected to contain:
#'   \itemize{
#'     \item \code{name}: rule name.
#'     \item \code{size_allowed}: integer size constraint.
#'     \item \code{match_fn}: function(\code{fam_df}) returning a scalar logical.
#'     \item \code{kept_rules}: list of character vectors (each element is an expression).
#'     \item \code{discarded_rules}: list of character vectors.
#'     \item \code{conflict_rules}: list of character vectors, each representing a candidate group.
#'     \item Optional: \code{alpha} (numeric), \code{kept_fn}, \code{discarded_fn}
#'           for graph-dependent deterministic selection.
#'   }
#'
#' @param seed Optional integer. If not \code{NULL}, sets the RNG seed before
#'   randomized choices.
#'
#' @param return_diagnostics Logical. If \code{FALSE} (default), returns only
#'   \code{kept} and \code{discarded}. If \code{TRUE}, returns additional fields
#'   describing intermediate sets and whether \code{igraph::largest_ivs()} was used.
#'
#' @param g_sub Optional \code{igraph} graph representing the same component as
#'   \code{fam_df}. Vertex names (\code{V(g_sub)$name}) must match \code{fam_df$.name}.
#'   Required to use the maximum independent set-based resolution (\code{rule$alpha > 1}).
#'
#' @details
#' Deterministic selection is computed via rule expressions evaluated over \code{fam_df}
#' columns (vectorized logical masks). Graph-dependent deterministic rules can be
#' provided via \code{rule$kept_fn} and \code{rule$discarded_fn}, which should return
#' character vectors of node names.
#'
#' Conflict groups are constructed from \code{rule$conflict_rules}, excluding any
#' nodes already classified a priori as kept or discarded.
#'
#' When \code{rule$alpha > 1} and \code{g_sub} is supplied, eligible nodes are
#' restricted to:
#' \itemize{
#'   \item all nodes excluding \code{discarded_apriori}, and
#'   \item optionally excluding \code{kept_apriori} and their 1-hop neighbors
#'         (to enforce compatibility with pre-fixed selections).
#' }
#' The function then samples uniformly among maximum independent vertex sets returned
#' by \code{igraph::largest_ivs()} on the induced subgraph.
#'
#' @return
#' If \code{return_diagnostics = FALSE}, a list with:
#' \itemize{
#'   \item \code{kept}: character vector of selected node names.
#'   \item \code{discarded}: character vector of discarded node names.
#' }
#'
#' If \code{return_diagnostics = TRUE}, a list with:
#' \itemize{
#'   \item \code{structure_name}
#'   \item \code{kept_apriori}, \code{discarded_apriori}
#'   \item \code{candidates}, \code{candidates_complex}
#'   \item \code{kept_afterrandomization}, \code{discarded_afterrandomization}
#'   \item \code{kept}, \code{discarded}
#'   \item \code{used_largest_ivs}
#' }
#'
#' @examples
#' \dontrun{
#' # fam_df: component-level node descriptors (size <= 5)
#' # rules_db: list of structural rules
#' out <- rkn_select_random_upto5(fam_df, rules_db, seed = 1)
#' out_diag <- rkn_select_random_upto5(fam_df, rules_db, seed = 1, return_diagnostics = TRUE, g_sub = g_sub)
#' }
#'
#' @export
rkn_select_random_upto5 <- function(fam_df,
                                    rules_db,
                                    seed = NULL,
                                    return_diagnostics = FALSE,
                                    g_sub = NULL) {
  if (!is.data.frame(fam_df)) {
    stop("`fam_df` must be a data.frame.")
  }
  if (!is.list(rules_db)) {
    stop("`rules_db` must be a list of rule objects.")
  }

  # Harmonize identifier columns.
  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name <- fam_df$.name
  if (is.null(fam_df$.name)) {
    stop("`fam_df` must contain a `.name` column (or `name` to be promoted).")
  }

  # Match the component to a rule (size <= 5).
  rule <- find_matching_rule(fam_df, rules_db)

  if (is.null(rule)) {
    kept_final <- character(0)
    discarded_final <- fam_df$.name

    if (!isTRUE(return_diagnostics)) {
      return(list(kept = kept_final, discarded = discarded_final))
    }

    return(list(
      structure_name = "unknown",
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidates = character(0),
      candidates_complex = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = FALSE
    ))
  }

  # Deterministic, table-driven selection.
  kept_apriori <- apply_rule_filter(fam_df, unlist(rule$kept_rules))
  discarded_apriori <- apply_rule_filter(fam_df, unlist(rule$discarded_rules))

  # Optional graph-dependent deterministic selection.
  if (!is.null(g_sub)) {
    if (!is.null(rule$kept_fn) && is.function(rule$kept_fn)) {
      kept_apriori <- unique(c(kept_apriori, rule$kept_fn(fam_df, g_sub)))
    }
    if (!is.null(rule$discarded_fn) && is.function(rule$discarded_fn)) {
      discarded_apriori <- unique(c(discarded_apriori, rule$discarded_fn(fam_df, g_sub)))
    }
  }

  # Build conflict groups excluding any a priori decisions.
  conflict_groups <- list()
  candidates_all <- character(0)

  if (length(rule$conflict_rules) > 0) {
    for (grp in rule$conflict_rules) {
      cand_nodes <- apply_rule_filter(fam_df, grp)
      cand_nodes <- setdiff(cand_nodes, union(kept_apriori, discarded_apriori))
      cand_nodes <- unique(cand_nodes)
      if (length(cand_nodes) > 0) {
        conflict_groups[[length(conflict_groups) + 1]] <- cand_nodes
      }
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
    # Eligible nodes: exclude discarded a priori and (optionally) neighbors of kept a priori.
    eligible <- setdiff(fam_df$.name, discarded_apriori)

    if (length(kept_apriori) > 0) {
      nb_idx <- igraph::neighborhood(g_sub, order = 1, nodes = kept_apriori)
      nb_names <- unique(setdiff(
        unlist(lapply(nb_idx, function(vs) igraph::V(g_sub)$name[vs])),
        kept_apriori
      ))
      eligible <- setdiff(eligible, union(kept_apriori, nb_names))
    }

    g_ind <- igraph::induced_subgraph(g_sub, vids = eligible)
    sets <- igraph::largest_ivs(g_ind)
    set_names <- lapply(sets, function(s) igraph::V(g_ind)$name[s])

    if (!is.null(seed)) set.seed(seed)
    kept_afterrandomization <- if (length(set_names) > 0) {
      set_names[[sample(seq_along(set_names), 1)]]
    } else {
      character(0)
    }

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

      kept_afterrandomization <- if (length(U) >= k) {
        sample(U, size = k, replace = FALSE)
      } else {
        chosen <- character(0)
        for (grp in conflict_groups) {
          avail <- setdiff(grp, chosen)
          if (length(avail) > 0) chosen <- c(chosen, sample(avail, 1))
        }
        unique(chosen)
      }
    }

    discarded_afterrandomization <- setdiff(candidates_all, kept_afterrandomization)
  }

  kept_afterrandomization <- unique(kept_afterrandomization)

  # Final reconciliation (kept overrides discarded).
  kept_final <- unique(c(kept_apriori, kept_afterrandomization))
  discarded_final <- unique(c(discarded_apriori, discarded_afterrandomization))
  discarded_final <- setdiff(discarded_final, kept_final)

  # Complexity flag for diagnostics (multiple groups or independent-set resolution).
  is_complex <- (length(conflict_groups) > 1) ||
    (!is.null(rule$alpha) && isTRUE(rule$alpha > 1))
  candidates_complex <- if (is_complex) unique(candidates_all) else character(0)

  if (!isTRUE(return_diagnostics)) {
    return(list(kept = kept_final, discarded = discarded_final))
  }

  if (!length(candidates_all)) candidates_all <- character(0)
  if (!length(discarded_afterrandomization)) discarded_afterrandomization <- character(0)
  if (!length(kept_afterrandomization)) kept_afterrandomization <- character(0)

  list(
    structure_name = rule$name,
    kept_apriori = unique(kept_apriori),
    discarded_apriori = unique(discarded_apriori),
    candidates = unique(candidates_all),
    candidates_complex = unique(candidates_complex),
    kept_afterrandomization = unique(kept_afterrandomization),
    discarded_afterrandomization = unique(discarded_afterrandomization),
    kept = kept_final,
    discarded = discarded_final,
    used_largest_ivs = used_largest_ivs
  )
}
