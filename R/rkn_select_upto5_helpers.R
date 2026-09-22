# Helpers shared between rkn_select_random_upto5() and rkn_select_criteria_upto5().
# The two selectors are identical up to the point where an ambiguous conflict group
# must be resolved (uniform random draw vs. criteria-based scoring); everything below
# is that shared, pre-resolution part.

#' Internal helper: build per-node diagnostics for upto5 selectors
#'
#' @description
#' Constructs the standard per-node diagnostics data frame (`pre_class`, `post_class`,
#' `is_selected`, `pre_reason`) shared by \code{rkn_select_random_upto5()} and
#' \code{rkn_select_criteria_upto5()}. If a node were to appear in more than one input
#' set (should not happen for well-formed calls, since \code{conflict_groups} already
#' excludes a priori decisions), the last of these wins, in order: \code{kept_apriori},
#' \code{discarded_apriori}, \code{candidate_simple}, \code{candidate_complex}.
#'
#' @keywords internal
#' @noRd
.rkn_build_node_diag_upto5 <- function(ids,
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

#' Internal helper: a priori decisions and conflict groups for upto5 selectors
#'
#' @description
#' Computes the deterministic (table-driven) a priori kept/discarded decisions from
#' \code{rule$kept_rules}/\code{rule$discarded_rules} (and, if \code{g_sub} is supplied,
#' the graph-dependent \code{rule$kept_fn}/\code{rule$discarded_fn}), then builds the
#' remaining ambiguous conflict groups from \code{rule$conflict_rules}, always excluding
#' nodes already decided a priori. Shared, resolution-agnostic step between
#' \code{rkn_select_random_upto5()} and \code{rkn_select_criteria_upto5()}.
#'
#' @return A list with \code{kept_apriori}, \code{discarded_apriori}, \code{conflict_groups}
#'   (a list of character vectors, one per group), and \code{candidates_all} (their union).
#'
#' @keywords internal
#' @noRd
.rkn_apriori_and_conflicts_upto5 <- function(fam_df, rule, g_sub = NULL) {
  kept_apriori      <- apply_rule_filter(fam_df, unlist(rule$kept_rules))
  discarded_apriori <- apply_rule_filter(fam_df, unlist(rule$discarded_rules))

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

  list(
    kept_apriori = kept_apriori,
    discarded_apriori = discarded_apriori,
    conflict_groups = conflict_groups,
    candidates_all = candidates_all
  )
}

#' Internal helper: exclude an already-kept set and its neighbours from eligibility
#'
#' @description
#' Removes \code{kept_apriori} and its 1-hop neighbours in \code{g_sub} from
#' \code{eligible}, so that a subsequent \code{igraph::largest_ivs()} call on the
#' induced subgraph of \code{eligible} cannot return a set adjacent to what is already
#' kept. Shared between \code{rkn_select_random_upto5()} and
#' \code{rkn_select_criteria_upto5()}'s `alpha`-based resolution branch.
#'
#' @keywords internal
#' @noRd
.rkn_exclude_neighbors_of_kept <- function(g_sub, eligible, kept_apriori) {
  if (length(kept_apriori) == 0) return(eligible)
  nb_idx  <- igraph::neighborhood(g_sub, order = 1, nodes = kept_apriori)
  nb_names <- unique(setdiff(
    unlist(lapply(nb_idx, function(vs) igraph::V(g_sub)$name[vs])),
    kept_apriori
  ))
  setdiff(eligible, union(kept_apriori, nb_names))
}
