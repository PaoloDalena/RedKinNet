# ---------- Internal constructor and shared helpers -------------------------

#' Internal constructor for `rkn_criterion` objects
#'
#' @keywords internal
#' @noRd
.rkn_new_criterion <- function(type, params, score_fn) {
  structure(
    list(type = type, params = params, score_fn = score_fn),
    class = c(paste0("rkn_criterion_", type), "rkn_criterion")
  )
}

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

#' Preference for the mean of a continuous covariate
#'
#' @description
#' Scores a candidate kept-set by how close the mean of `covariate`, computed over
#' the nodes in the kept-set, is to a target value. The score is 1 when the mean
#' equals `target` and decreases linearly to 0 at a distance of `tolerance`.
#'
#' @param covariate Character. Name of a numeric column expected in `fam_df`.
#' @param target Numeric. Target value for the mean of `covariate` within the kept-set.
#' @param tolerance Numeric (> 0). Distance from `target` at which the score reaches 0.
#'   There is no default: the appropriate scale depends on the covariate and must be
#'   supplied explicitly, since within-family variability is not a reliable estimator
#'   of scale for the small family sizes handled by \code{rkn_select_criteria_upto5()}.
#'
#' @return An object of class \code{c("rkn_criterion_target_mean", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_criterion_cutoff}}, \code{\link{rkn_criterion_level}},
#'   \code{\link{rkn_criterion_category_balance}}
#' @export
rkn_criterion_target_mean <- function(covariate, target, tolerance) {
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.")
  }
  if (!is.numeric(target) || length(target) != 1L || is.na(target)) {
    stop("`target` must be a single non-missing numeric value.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be a single positive numeric value.")
  }

  score_fn <- function(fam_df, kept_ids) {
    x <- .rkn_get_covariate(fam_df, kept_ids, covariate)
    if (length(x) == 0) return(0)
    m <- mean(x, na.rm = TRUE)
    if (is.nan(m)) return(0)
    max(0, 1 - abs(m - target) / tolerance)
  }

  .rkn_new_criterion(
    type = "target_mean",
    params = list(covariate = covariate, target = target, tolerance = tolerance),
    score_fn = score_fn
  )
}

#' Preference for values above or below a cutoff on a continuous covariate
#'
#' @description
#' Scores a candidate kept-set by the proportion of its nodes whose `covariate` value
#' satisfies the preferred direction relative to `cutoff`.
#'
#' @param covariate Character. Name of a numeric column expected in `fam_df`.
#' @param cutoff Numeric. Threshold value.
#' @param prefer Character. Either `"above"` (prefer \code{covariate > cutoff}) or
#'   `"below"` (prefer \code{covariate < cutoff}).
#'
#' @return An object of class \code{c("rkn_criterion_cutoff", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_criterion_target_mean}}, \code{\link{rkn_criterion_level}},
#'   \code{\link{rkn_criterion_category_balance}}
#' @export
rkn_criterion_cutoff <- function(covariate, cutoff, prefer = c("above", "below")) {
  prefer <- match.arg(prefer)
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a single non-missing numeric value.")
  }

  score_fn <- function(fam_df, kept_ids) {
    x <- .rkn_get_covariate(fam_df, kept_ids, covariate)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(0)
    satisfies <- if (identical(prefer, "above")) x > cutoff else x < cutoff
    mean(satisfies)
  }

  .rkn_new_criterion(
    type = "cutoff",
    params = list(covariate = covariate, cutoff = cutoff, prefer = prefer),
    score_fn = score_fn
  )
}

#' Preference for a specific level of a covariate
#'
#' @description
#' Scores a candidate kept-set by the proportion of its nodes whose `covariate` value
#' equals `level`. Intended for dichotomous covariates, but works for any covariate
#' type supporting `==`.
#'
#' @param covariate Character. Name of a column expected in `fam_df`.
#' @param level The preferred value of `covariate` (matched with `==`).
#'
#' @return An object of class \code{c("rkn_criterion_level", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_criterion_target_mean}}, \code{\link{rkn_criterion_cutoff}},
#'   \code{\link{rkn_criterion_category_balance}}
#' @export
rkn_criterion_level <- function(covariate, level) {
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.")
  }
  if (length(level) != 1L || is.na(level)) {
    stop("`level` must be a single non-missing value.")
  }

  score_fn <- function(fam_df, kept_ids) {
    x <- .rkn_get_covariate(fam_df, kept_ids, covariate)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(0)
    mean(x == level)
  }

  .rkn_new_criterion(
    type = "level",
    params = list(covariate = covariate, level = level),
    score_fn = score_fn
  )
}

#' Preference for preserving categorical balance
#'
#' @description
#' Scores a candidate kept-set by 1 minus the total variation distance between the
#' category proportions of `covariate` within the kept-set and a reference
#' distribution. Total variation distance is bounded in `[0, 1]`, so the resulting
#' score is as well.
#'
#' @param covariate Character. Name of a column expected in `fam_df`.
#' @param target_props Optional named numeric vector of reference proportions
#'   (must sum to 1), one entry per targeted category level, e.g.
#'   \code{c(male = 0.5, female = 0.5)}. Categories absent from `target_props` are not
#'   assigned an implicit target of zero: the comparison is restricted to the named
#'   categories only, using their absolute proportion within the kept-set (so kept-set
#'   members falling in an unlisted category still dilute these proportions, but their
#'   mere presence is not itself penalized against an explicit zero target).
#'   If `NULL` (default), the reference distribution is the family's own observed
#'   category proportions in `fam_df` (i.e. "preserve this family's natural mix"); in
#'   that case every observed level already has a real reference proportion, so there
#'   are no omitted categories to begin with.
#'
#' @return An object of class \code{c("rkn_criterion_category_balance", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_criterion_target_mean}}, \code{\link{rkn_criterion_cutoff}},
#'   \code{\link{rkn_criterion_level}}
#' @export
rkn_criterion_category_balance <- function(covariate, target_props = NULL) {
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.")
  }
  if (!is.null(target_props)) {
    if (!is.numeric(target_props) || is.null(names(target_props)) || any(names(target_props) == "")) {
      stop("`target_props` must be a named numeric vector.")
    }
    if (abs(sum(target_props) - 1) > 1e-8) {
      stop("`target_props` must sum to 1.")
    }
  }

  score_fn <- function(fam_df, kept_ids) {
    x <- .rkn_get_covariate(fam_df, kept_ids, covariate)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(0)

    kept_props <- table(x) / length(x)

    if (is.null(target_props)) {
      all_x <- fam_df[[covariate]]
      all_x <- all_x[!is.na(all_x)]
      ref_props <- table(all_x) / length(all_x)

      lvls <- union(names(kept_props), names(ref_props))
      kp <- setNames(rep(0, length(lvls)), lvls)
      rp <- setNames(rep(0, length(lvls)), lvls)
      kp[names(kept_props)] <- as.numeric(kept_props)
      rp[names(ref_props)] <- as.numeric(ref_props)
    } else {
      # Compare only the explicitly targeted categories, using their absolute
      # proportion within the kept-set. A category missing from target_props is not
      # assigned an implicit target of zero: it simply is not part of the comparison,
      # though its members still dilute the proportions of the targeted categories.
      lvls <- names(target_props)
      kp <- setNames(rep(0, length(lvls)), lvls)
      matched <- intersect(names(kept_props), lvls)
      kp[matched] <- as.numeric(kept_props[matched])
      rp <- as.numeric(target_props)
    }

    tv_dist <- 0.5 * sum(abs(kp - rp))
    1 - tv_dist
  }

  .rkn_new_criterion(
    type = "category_balance",
    params = list(covariate = covariate, target_props = target_props),
    score_fn = score_fn
  )
}

#' @export
print.rkn_criterion <- function(x, ...) {
  cat(sprintf("<rkn_criterion: %s>\n", x$type))
  for (nm in names(x$params)) {
    val <- x$params[[nm]]
    val_str <- if (is.null(val)) "NULL" else paste(format(val), collapse = ", ")
    cat(sprintf("  %s: %s\n", nm, val_str))
  }
  invisible(x)
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
