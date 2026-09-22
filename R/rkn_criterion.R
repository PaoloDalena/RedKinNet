# ---------- Internal constructor -------------------------------------------

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
