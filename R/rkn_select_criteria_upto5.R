#' Criteria-based selection within connected components of size <= 5
#'
#' @description
#' Applies the RedKinNet node selection strategy to a connected component with
#' at most 5 nodes, resolving ambiguous conflict groups by score under one or
#' two covariate-based criteria instead of uniformly at random.
#'
#' Rule matching and deterministic (a priori) keep/discard decisions are
#' identical to \code{\link{rkn_select_random_upto5}} (same \code{rules_db}/
#' \code{custom_rules_db}, same \code{kept_rules}/\code{discarded_rules}/
#' \code{kept_fn}/\code{discarded_fn}). The two functions differ only in how
#' the remaining ambiguous conflict group is resolved:
#' \itemize{
#'   \item If \code{rule$alpha > 1} and \code{g_sub} is provided, every maximum
#'   independent set returned by \code{igraph::largest_ivs()} on the induced
#'   subgraph of eligible vertices is treated as a candidate, and the
#'   best-scoring one is kept (see \code{\link{.rkn_select_among_candidates}}).
#'   \item Otherwise, each single node in the (unique) conflict group is
#'   treated as a candidate, and the best-scoring one is kept.
#' }
#' Exact ties in the criterion score(s) are broken with a seeded random draw
#' restricted to the tied candidates.
#'
#' Rules with more than one disjoint conflict group (relevant only for a
#' \code{custom_rules_db}; no rule in the packaged \code{rules_db} has more
#' than one) are not supported and raise an error, since scoring their joint
#' resolution correctly requires enumerating combinations across groups, which
#' is not yet implemented.
#'
#' @inheritParams rkn_select_random_upto5
#' @param criterion_primary An \code{rkn_criterion} object (see
#'   \code{\link{rkn_criterion_target_mean}} and siblings) used to score
#'   candidate resolutions of the conflict group.
#' @param criterion_secondary Optional \code{rkn_criterion} object used to
#'   break ties in \code{criterion_primary}.
#' @param seed Optional integer RNG seed, used only to break exact ties in the
#'   criterion score(s).
#'
#' @return Same structure as \code{\link{rkn_select_random_upto5}}.
#'
#' @seealso \code{\link{rkn_select_random_upto5}}, \code{\link{rkn_criterion_target_mean}},
#'   \code{\link{rkn_score_candidate_set}}
#'
#' @export
rkn_select_criteria_upto5 <- function(fam_df,
                                      criterion_primary,
                                      criterion_secondary = NULL,
                                      custom_rules_db = NULL,
                                      seed = NULL,
                                      return_diagnostics = FALSE,
                                      g_sub = NULL) {

  if (missing(criterion_primary) || !inherits(criterion_primary, "rkn_criterion")) {
    stop("`criterion_primary` must be an rkn_criterion object (see rkn_criterion_*()).")
  }
  if (!is.null(criterion_secondary) && !inherits(criterion_secondary, "rkn_criterion")) {
    stop("`criterion_secondary` must be an rkn_criterion object or NULL.")
  }

  # Harmonize identifier columns.
  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  # Choose the rule database: internal `rules_db` (default) or user-supplied `custom_rules_db`.
  if (is.null(custom_rules_db)) {
    if (!exists("rules_db", inherits = TRUE)) {
      stop("Internal `rules_db` was not found. Reinstall/reload the package, or supply `custom_rules_db`.")
    }
    rules_db_use <- rules_db
  } else {
    rules_db_use <- custom_rules_db
  }

  # Match the component to a rule.
  rule <- find_matching_rule(fam_df, rules_db_use)

  # Build per-node diagnostics (aligned with rkn_select_random_upto5's fields).
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
      kept_aftercriteria = character(0),
      discarded_aftercriteria = character(0),
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

  if (length(conflict_groups) > 1) {
    stop(
      "rkn_select_criteria_upto5() does not support rules with more than one ",
      "disjoint conflict group (matched rule: '", rule$name, "'). ",
      "Scoring their joint resolution requires enumerating combinations across ",
      "groups, which is not yet implemented."
    )
  }

  # Criteria-based resolution.
  used_largest_ivs <- FALSE
  kept_aftercriteria <- character(0)
  discarded_aftercriteria <- character(0)

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
    candidate_sets <- lapply(sets, function(s) igraph::V(g_ind)$name[s])

    if (length(candidates_all) == 0) candidates_all <- eligible

    if (length(candidate_sets) > 0) {
      picked <- .rkn_select_among_candidates(
        candidate_sets = candidate_sets,
        fam_df = fam_df,
        fixed_kept_ids = kept_apriori,
        criterion_primary = criterion_primary,
        criterion_secondary = criterion_secondary,
        seed = seed
      )
      kept_aftercriteria <- unique(picked$chosen_set)
    }

    discarded_aftercriteria <- setdiff(candidates_all, kept_aftercriteria)
    used_largest_ivs <- TRUE

  } else if (length(conflict_groups) == 1) {

    grp <- conflict_groups[[1]]
    candidate_sets <- as.list(grp)

    picked <- .rkn_select_among_candidates(
      candidate_sets = candidate_sets,
      fam_df = fam_df,
      fixed_kept_ids = kept_apriori,
      criterion_primary = criterion_primary,
      criterion_secondary = criterion_secondary,
      seed = seed
    )
    kept_aftercriteria <- unique(picked$chosen_set)
    discarded_aftercriteria <- setdiff(candidates_all, kept_aftercriteria)
  }

  kept_aftercriteria <- unique(kept_aftercriteria)

  # Final reconciliation (kept overrides discarded).
  kept_final <- unique(c(kept_apriori, kept_aftercriteria))
  discarded_final <- unique(c(discarded_apriori, discarded_aftercriteria))
  discarded_final <- setdiff(discarded_final, kept_final)

  # Candidate complexity flag (complex if resolved via largest_ivs).
  is_complex <- !is.null(rule$alpha) && isTRUE(rule$alpha > 1)

  candidate_complex <- if (is_complex) unique(candidates_all) else character(0)
  candidate_simple  <- if (!is_complex) unique(candidates_all) else character(0)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  }

  # Minimal pre-reason map for diagnostics/plotting.
  pre_reason_map <- character(0)
  if (used_largest_ivs && length(candidate_complex) > 0) {
    pre_reason_map[candidate_complex] <- "criteria_largest_ivs"
  } else if (length(conflict_groups) == 1 && length(candidate_simple) > 0) {
    pre_reason_map[candidate_simple] <- "criteria_single_conflict_group"
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
    kept_aftercriteria = kept_aftercriteria,
    discarded_aftercriteria = unique(discarded_aftercriteria),
    kept = kept_final,
    discarded = discarded_final,
    used_largest_ivs = used_largest_ivs,
    node_diagnostics = node_diagnostics,
    family_diagnostics = family_diagnostics
  )
}
