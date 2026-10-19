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
#' Scores a candidate set of individuals by how close the mean of `covariate` over its
#' members is to a target value. The score is 1 when the mean equals `target` and
#' decreases linearly to 0 at a distance of `tolerance`.
#'
#' Criteria are used by \code{\link{rkn_select_criteria}}: within each family, where
#' the structural rules leave several sets of unrelated individuals to choose from,
#' each candidate set is scored and the best one is selected.
#'
#' @param covariate Character. Name of a numeric vertex attribute of the network, such
#'   as \code{"age"} for \code{V(g)$age}.
#' @param target Numeric. Target value for the mean of `covariate` within the candidate
#'   set.
#' @param tolerance Numeric (> 0). Distance from `target` at which the score reaches 0.
#'   There is no default: the appropriate scale depends on the covariate and must be
#'   supplied explicitly, since the variability within a family, which has only a few
#'   members, is not a reliable estimate of scale.
#'
#' @return An object of class \code{c("rkn_criterion_target_mean", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_select_criteria}}, \code{\link{rkn_criterion_cutoff}},
#'   \code{\link{rkn_criterion_level}}, \code{\link{rkn_criterion_category_balance}}
#'
#' @examples
#' rkn_criterion_target_mean("age", target = 40, tolerance = 20)
#'
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
#' Scores a candidate set of individuals by the proportion of its members whose
#' `covariate` value lies on the preferred side of `cutoff`.
#'
#' Criteria are used by \code{\link{rkn_select_criteria}}: within each family, where
#' the structural rules leave several sets of unrelated individuals to choose from,
#' each candidate set is scored and the best one is selected.
#'
#' @param covariate Character. Name of a numeric vertex attribute of the network, such
#'   as \code{"age"} for \code{V(g)$age}.
#' @param cutoff Numeric. Threshold value.
#' @param prefer Character. Either `"above"` (prefer \code{covariate > cutoff}) or
#'   `"below"` (prefer \code{covariate < cutoff}).
#'
#' @return An object of class \code{c("rkn_criterion_cutoff", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_select_criteria}}, \code{\link{rkn_criterion_target_mean}},
#'   \code{\link{rkn_criterion_level}}, \code{\link{rkn_criterion_category_balance}}
#'
#' @examples
#' # Prefer individuals younger than 50
#' rkn_criterion_cutoff("age", cutoff = 50, prefer = "below")
#'
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
#' Scores a candidate set of individuals by the proportion of its members whose
#' `covariate` value equals `level`. Intended for dichotomous covariates, but works for
#' any covariate type supporting `==`.
#'
#' Criteria are used by \code{\link{rkn_select_criteria}}: within each family, where
#' the structural rules leave several sets of unrelated individuals to choose from,
#' each candidate set is scored and the best one is selected.
#'
#' @param covariate Character. Name of a vertex attribute of the network, such as
#'   \code{"sex"} for \code{V(g)$sex}.
#' @param level The preferred value of `covariate` (matched with `==`).
#'
#' @return An object of class \code{c("rkn_criterion_level", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_select_criteria}}, \code{\link{rkn_criterion_target_mean}},
#'   \code{\link{rkn_criterion_cutoff}}, \code{\link{rkn_criterion_category_balance}}
#'
#' @examples
#' # Prefer women
#' rkn_criterion_level("sex", "F")
#'
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
#' Scores a candidate set of individuals by 1 minus the total variation distance between
#' the category proportions of `covariate` among its members and a reference
#' distribution. Total variation distance lies between 0 and 1, and so does the score.
#'
#' Criteria are used by \code{\link{rkn_select_criteria}}: within each family, where
#' the structural rules leave several sets of unrelated individuals to choose from,
#' each candidate set is scored and the best one is selected.
#'
#' @param covariate Character. Name of a vertex attribute of the network, such as
#'   \code{"sex"} for \code{V(g)$sex}.
#' @param target_props Optional named numeric vector of reference proportions
#'   (must sum to 1), one entry per targeted category level, e.g.
#'   \code{c(male = 0.5, female = 0.5)}. Categories absent from `target_props` are not
#'   assigned an implicit target of zero: the comparison is restricted to the named
#'   categories only, using their absolute proportion within the candidate set (so
#'   members falling in an unlisted category still dilute these proportions, but their
#'   mere presence is not itself penalized against an explicit zero target).
#'   If `NULL` (default), the reference distribution is the observed category
#'   proportions of the whole family (i.e. "preserve this family's natural mix"); in
#'   that case every observed level already has a real reference proportion, so there
#'   are no omitted categories to begin with.
#'
#' @return An object of class \code{c("rkn_criterion_category_balance", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_select_criteria}}, \code{\link{rkn_criterion_target_mean}},
#'   \code{\link{rkn_criterion_cutoff}}, \code{\link{rkn_criterion_level}}
#'
#' @examples
#' # Keep the mix of each family
#' rkn_criterion_category_balance("sex")
#'
#' # Aim at equal proportions of women and men
#' rkn_criterion_category_balance("sex", c(F = 0.5, M = 0.5))
#'
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
      kp <- stats::setNames(rep(0, length(lvls)), lvls)
      rp <- stats::setNames(rep(0, length(lvls)), lvls)
      kp[names(kept_props)] <- as.numeric(kept_props)
      rp[names(ref_props)] <- as.numeric(ref_props)
    } else {
      # Compare only the explicitly targeted categories, using their absolute
      # proportion within the kept-set. A category missing from target_props is not
      # assigned an implicit target of zero: it simply is not part of the comparison,
      # though its members still dilute the proportions of the targeted categories.
      lvls <- names(target_props)
      kp <- stats::setNames(rep(0, length(lvls)), lvls)
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
