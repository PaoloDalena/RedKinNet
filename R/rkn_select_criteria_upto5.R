#' Criteria-based selection within connected components of size <= 5
#'
#' @description
#' Applies the RedKinNet node selection strategy to a connected component with
#' at most 5 nodes, resolving ambiguous conflict groups by score under one or
#' two covariate-based criteria instead of uniformly at random.
#'
#' Rule matching and deterministic (a priori) keep/discard decisions are
#' identical to \code{.rkn_select_random_upto5()} (same \code{rules_db}/
#' \code{custom_rules_db}, same \code{kept_rules}/\code{discarded_rules}/
#' \code{kept_fn}/\code{discarded_fn}). The two functions differ only in how
#' the remaining ambiguous conflict group is resolved:
#' \itemize{
#'   \item If \code{rule$use_largest_ivs} is \code{TRUE} and \code{g_sub} is
#'   provided, every maximum independent set returned by \code{igraph::largest_ivs()}
#'   on the induced subgraph of eligible vertices is treated as a candidate, and
#'   the candidate with the best primary score is kept (ties are broken by the
#'   secondary criterion, if supplied).
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
#' @inheritParams .rkn_select_random_upto5
#' @param criterion_primary An \code{rkn_criterion} object (see
#'   \code{rkn_criterion_level()} and \code{rkn_criterion_cutoff()}) used to score
#'   candidate resolutions of the conflict group.
#' @param criterion_secondary Optional \code{rkn_criterion} object used to
#'   break ties in \code{criterion_primary}.
#' @param seed Optional integer RNG seed, used only to break exact ties in the
#'   criterion score(s).
#'
#' @return Same structure as \code{.rkn_select_random_upto5()}.
#'
#' @seealso \code{.rkn_select_random_upto5()}, \code{rkn_criterion_level()},
#'   \code{.rkn_score_candidate_set()}
#'
#' @keywords internal
#' @noRd
.rkn_select_criteria_upto5 <- function(fam_df,
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

    node_diagnostics <- .rkn_df(
      id = fam_df$.name,
      pre_class = rep("discarded_apriori", nrow(fam_df)),
      post_class = rep("discarded", nrow(fam_df)),
      is_selected = rep(FALSE, nrow(fam_df)),
      pre_reason = rep("no_matching_rule", nrow(fam_df))
    )

    family_diagnostics <- .rkn_df(
      n_nodes = nrow(fam_df),
      mis_method = NA_character_,
      structure_name = "unknown"
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
      mis_method = NA_character_,
      node_diagnostics = node_diagnostics,
      family_diagnostics = family_diagnostics
    ))
  }

  ac <- .rkn_apriori_and_conflicts_upto5(fam_df, rule, g_sub)
  kept_apriori      <- ac$kept_apriori
  discarded_apriori <- ac$discarded_apriori
  conflict_groups   <- ac$conflict_groups
  candidates_all    <- ac$candidates_all

  if (length(conflict_groups) > 1) {
    stop(
      ".rkn_select_criteria_upto5() does not support rules with more than one ",
      "disjoint conflict group (matched rule: '", rule$name, "'). ",
      "Scoring their joint resolution requires enumerating combinations across ",
      "groups, which is not yet implemented."
    )
  }

  used_largest_ivs <- FALSE
  kept_aftercriteria <- character(0)
  discarded_aftercriteria <- character(0)

  if (isTRUE(rule$use_largest_ivs) && !is.null(g_sub)) {

    # Nodes adjacent to an a-priori kept node can no longer be selected.
    eligible <- setdiff(fam_df$.name, discarded_apriori)
    eligible <- .rkn_exclude_neighbors_of_kept(g_sub, eligible, kept_apriori)

    g_ind <- igraph::induced_subgraph(g_sub, vids = eligible)
    sets  <- igraph::largest_ivs(g_ind)
    ind_names <- igraph::V(g_ind)$name
    candidate_sets <- vector("list", length(sets))
    for (j in seq_along(sets)) candidate_sets[[j]] <- ind_names[sets[[j]]]

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

  # Kept overrides discarded.
  kept_final <- unique(c(kept_apriori, kept_aftercriteria))
  discarded_final <- unique(c(discarded_apriori, discarded_aftercriteria))
  discarded_final <- setdiff(discarded_final, kept_final)

  is_complex <- isTRUE(rule$use_largest_ivs)

  candidate_complex <- if (is_complex) unique(candidates_all) else character(0)
  candidate_simple  <- if (!is_complex) unique(candidates_all) else character(0)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  }

  pre_reason_map <- character(0)
  if (used_largest_ivs && length(candidate_complex) > 0) {
    pre_reason_map[candidate_complex] <- "criteria_largest_ivs"
  } else if (length(conflict_groups) == 1 && length(candidate_simple) > 0) {
    pre_reason_map[candidate_simple] <- "criteria_single_conflict_group"
  }

  if (used_largest_ivs) {
    split <- .rkn_split_enumerated_candidates(candidate_complex, candidate_sets)
    kept_apriori <- c(kept_apriori, split$kept)
    discarded_apriori <- c(discarded_apriori, split$discarded)
    candidate_complex <- split$complex
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

  family_diagnostics <- .rkn_df(
    n_nodes = nrow(fam_df),
    mis_method = if (used_largest_ivs) "enumeration" else NA_character_,
    structure_name = rule$name
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
    mis_method = if (used_largest_ivs) "enumeration" else NA_character_,
    node_diagnostics = node_diagnostics,
    family_diagnostics = family_diagnostics
  )
}
