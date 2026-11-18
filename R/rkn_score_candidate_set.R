#' Internal helper: extract covariate values for a candidate kept-set
#'
#' @keywords internal
#' @noRd
.rkn_get_covariate <- function(fam_df, kept_ids, covariate) {
  if (is.null(fam_df[[covariate]])) {
    stop(sprintf("Covariate \"%s\" not found.", covariate))
  }
  idx <- fam_df$.name %in% kept_ids
  fam_df[[covariate]][idx]
}

#' Score a candidate kept-set against one or two covariate criteria
#'
#' @description
#' Evaluates `criterion_primary` (and, if supplied, `criterion_secondary`) against a
#' fully-specified candidate kept-set. Intended both as the scoring step used
#' internally by \code{rkn_select_criteria_upto5()}/\code{rkn_select_criteria_6plus()},
#' and as a standalone diagnostic tool for reporting how well an already-decided
#' selection satisfies a set of criteria.
#'
#' @param fam_df A data frame describing the nodes in a family. Must contain a `.name`
#'   column and the covariate column(s) referenced by the supplied criteria.
#' @param kept_ids Character vector of node identifiers in the candidate kept-set.
#' @param criterion_primary An `rkn_criterion` object (see \code{rkn_criterion_level()} and
#'   \code{rkn_criterion_cutoff()}).
#' @param criterion_secondary Optional `rkn_criterion` object.
#'
#' @return A one-row data frame with columns `primary_type`, `primary_score`,
#'   `secondary_type`, `secondary_score` (the latter two `NA` if `criterion_secondary`
#'   is not supplied).
#'
#' @keywords internal
#' @noRd
rkn_score_candidate_set <- function(fam_df, kept_ids, criterion_primary, criterion_secondary = NULL) {
  if (missing(criterion_primary) || !inherits(criterion_primary, "rkn_criterion")) {
    stop("`criterion_primary` must be an rkn_criterion object (see rkn_criterion_*()).")
  }
  if (!is.null(criterion_secondary) && !inherits(criterion_secondary, "rkn_criterion")) {
    stop("`criterion_secondary` must be an rkn_criterion object or NULL.")
  }

  primary_score <- criterion_primary$score_fn(fam_df, kept_ids)
  secondary_score <- if (is.null(criterion_secondary)) {
    NA_real_
  } else {
    criterion_secondary$score_fn(fam_df, kept_ids)
  }

  data.frame(
    primary_type = criterion_primary$type,
    primary_score = primary_score,
    secondary_type = if (is.null(criterion_secondary)) NA_character_ else criterion_secondary$type,
    secondary_score = secondary_score,
    stringsAsFactors = FALSE
  )
}

#' Internal helper: pick one candidate node-set using up to two criteria
#'
#' @description
#' Scores every candidate in `candidate_sets` against `criterion_primary` (evaluated
#' over the union of `fixed_kept_ids` and each candidate's own contribution), keeps
#' only the candidates achieving the best primary score, and -- if more than one
#' remains and `criterion_secondary` is supplied -- narrows further using the
#' secondary score. Any remaining exact tie is broken with a seeded random draw
#' restricted to the tied candidates.
#'
#' @keywords internal
#' @noRd
.rkn_select_among_candidates <- function(candidate_sets,
                                         fam_df,
                                         fixed_kept_ids = character(0),
                                         criterion_primary,
                                         criterion_secondary = NULL,
                                         seed = NULL) {
  if (!is.null(seed)) withr::local_preserve_seed()

  n_cand <- length(candidate_sets)
  if (n_cand == 0) stop("`candidate_sets` must be a non-empty list.")
  if (n_cand == 1) return(list(chosen_index = 1L, chosen_set = candidate_sets[[1]]))

  primary_scores <- numeric(n_cand)
  for (j in seq_len(n_cand)) {
    full_kept <- unique(c(fixed_kept_ids, candidate_sets[[j]]))
    primary_scores[j] <- criterion_primary$score_fn(fam_df, full_kept)
  }

  best_primary <- max(primary_scores)
  tied_idx <- which(primary_scores == best_primary)

  if (length(tied_idx) > 1 && !is.null(criterion_secondary)) {
    secondary_scores <- numeric(length(tied_idx))
    for (k in seq_along(tied_idx)) {
      full_kept <- unique(c(fixed_kept_ids, candidate_sets[[tied_idx[k]]]))
      secondary_scores[k] <- criterion_secondary$score_fn(fam_df, full_kept)
    }
    best_secondary <- max(secondary_scores)
    tied_idx <- tied_idx[secondary_scores == best_secondary]
  }

  if (length(tied_idx) > 1) {
    if (!is.null(seed)) set.seed(seed)
    chosen <- tied_idx[sample(length(tied_idx), 1)]
  } else {
    chosen <- tied_idx[1]
  }

  list(chosen_index = chosen, chosen_set = candidate_sets[[chosen]])
}

#' Internal helper: per-individual contributions of an additive criterion
#'
#' The \code{level} and \code{cutoff} criteria score a set by the proportion of its
#' members, among those with a non-missing covariate, that meet a condition. Among sets
#' of the same size whose members all have the covariate, this proportion is ranked as
#' the number of members meeting the condition, a sum over individuals that the exact
#' engine can maximise. Returns that 0/1 contribution for each of \code{ids}, or
#' \code{NULL} for the other criteria or when some of \code{ids} lack the covariate.
#'
#' @keywords internal
#' @noRd
.rkn_additive_contribution <- function(criterion, fam_df, ids) {
  if (!criterion$type %in% c("level", "cutoff")) return(NULL)
  x <- fam_df[[criterion$params$covariate]]
  if (is.null(x)) return(NULL)
  x <- x[match(ids, fam_df$.name)]
  if (anyNA(x)) return(NULL)
  hit <- if (identical(criterion$type, "level")) {
    x == criterion$params$level
  } else if (identical(criterion$params$prefer, "above")) {
    x > criterion$params$cutoff
  } else {
    x < criterion$params$cutoff
  }
  as.numeric(hit)
}

#' Internal helper: choose a maximum independent set of a graph by criteria
#'
#' Additive criteria (see \code{.rkn_additive_contribution()}) are optimised exactly by
#' the exact engine, with integer vertex weights that rank sets by size, then by the
#' primary and then by the secondary criterion; the remaining ties are broken uniformly
#' at random. A criterion that is not additive is evaluated on candidate sets: all the
#' maximum sets (optimal for the preceding additive criterion, if any) when there are at
#' most \code{limit}, otherwise a uniform sample of \code{limit} of them, which is
#' reported as \code{"exact_sample"}. A graph too wide for the engine is first reduced by
#' dominance (\code{"exact_reduced"}), as in \code{.rkn_draw_maximum_set()}; the
#' reduction uses the weights, so the additive criteria remain optimised exactly, while
#' the others are evaluated on the maximum sets of the reduced graph.
#'
#' @param g An igraph object with vertex names: the individuals to resolve.
#' @param fixed_kept Individuals of the family already kept, included in the scores.
#' @param error_if_wide If \code{FALSE}, a graph still too wide after the reduction
#'   returns \code{NULL} instead of an error.
#' @return A list with \code{kept}, the names of the selected vertices; \code{method};
#'   and \code{searched}, the names of the vertices of the graph searched.
#'
#' @keywords internal
#' @noRd
.rkn_criteria_maximum_set <- function(g, fam_df, fixed_kept, criterion_primary, criterion_secondary,
                                      seed = NULL, max_width = 20, limit = 10000, error_if_wide = TRUE) {
  if (!is.null(seed)) {
    withr::local_preserve_seed()
    set.seed(seed)
  }
  ids <- igraph::V(g)$name
  A <- .rkn_adjacency(g)
  base <- length(ids) + 1
  c1 <- .rkn_additive_contribution(criterion_primary, fam_df, ids)
  c2 <- if (!is.null(c1) && !is.null(criterion_secondary)) {
    .rkn_additive_contribution(criterion_secondary, fam_df, ids)
  }
  fully_additive <- !is.null(c1) && (is.null(criterion_secondary) || !is.null(c2))
  weight <- if (is.null(c1)) {
    rep(1, length(ids))
  } else if (fully_additive && !is.null(c2)) {
    base^2 + base * c1 + c2
  } else {
    base + c1
  }

  keep <- rep(TRUE, length(ids))
  method <- "exact"
  sol <- .rkn_mis_exact(A, max_width, weight)
  if (!sol$feasible) {
    keep <- .rkn_dominance_reduce(A, weight)
    sol <- .rkn_mis_exact(A[keep, keep, drop = FALSE], max_width, weight[keep])
    method <- "exact_reduced"
  }
  if (!sol$feasible) {
    if (!error_if_wide) return(NULL)
    stop(.rkn_intricate_message(ids, sol$width, max_width), call. = FALSE)
  }
  ids <- ids[keep]
  if (fully_additive) return(list(kept = ids[.rkn_mis_sample(sol)], method = method, searched = ids))

  sets <- .rkn_mis_enumerate(sol, limit)
  if (is.null(sets)) {
    sets <- .rkn_mis_sample_many(sol, limit)
    sets <- sets[, !duplicated(t(sets)), drop = FALSE]
    if (identical(method, "exact")) method <- "exact_sample"
  }
  candidate_sets <- vector("list", ncol(sets))
  for (j in seq_len(ncol(sets))) candidate_sets[[j]] <- ids[sets[, j]]
  picked <- .rkn_select_among_candidates(
    candidate_sets = candidate_sets, fam_df = fam_df, fixed_kept_ids = fixed_kept,
    criterion_primary = criterion_primary, criterion_secondary = criterion_secondary, seed = seed
  )
  list(kept = picked$chosen_set, method = method, searched = ids)
}
