# ---------- Internal helper --------------------------------------------------

#' Internal helper: extract covariate values for a candidate kept-set
#'
#' @keywords internal
#' @noRd
.rkn_get_covariate <- function(fam_df, kept_ids, covariate) {
  if (is.null(fam_df[[covariate]])) {
    stop(sprintf("Covariate '%s' not found in fam_df.", covariate))
  }
  idx <- fam_df$.name %in% kept_ids
  fam_df[[covariate]][idx]
}

# ---------- Scoring engine ---------------------------------------------------

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
#' @param criterion_primary An `rkn_criterion` object (see \code{\link{rkn_criterion_target_mean}}
#'   and siblings).
#' @param criterion_secondary Optional `rkn_criterion` object.
#'
#' @return A one-row data frame with columns `primary_type`, `primary_score`,
#'   `secondary_type`, `secondary_score` (the latter two `NA` if `criterion_secondary`
#'   is not supplied).
#'
#' @export
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
