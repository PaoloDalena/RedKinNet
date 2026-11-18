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

#' Prefer individuals above or below a cutoff
#'
#' @description
#' Scores a candidate set of individuals by the proportion of its members whose
#' `covariate` value lies on the preferred side of `cutoff`.
#'
#' Criteria are used by \code{\link{rkn_select_criteria}}: within each family, where
#' several sets of unrelated individuals remain to choose from, each candidate set is
#' scored and the best one is selected.
#'
#' @param covariate Character. Name of a numeric vertex attribute of the network, such
#'   as \code{"age"} for \code{V(g)$age}.
#' @param cutoff Numeric. Threshold value.
#' @param prefer Character. Either `"above"` (prefer \code{covariate > cutoff}) or
#'   `"below"` (prefer \code{covariate < cutoff}).
#'
#' @return An object of class \code{c("rkn_criterion_cutoff", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_select_criteria}}, \code{\link{rkn_criterion_level}};
#'   \code{\link{rkn_select_global}} for a target composition of the whole selection.
#'
#' @examples
#' # Prefer individuals younger than 50
#' rkn_criterion_cutoff("age", cutoff = 50, prefer = "below")
#'
#' @export
rkn_criterion_cutoff <- function(covariate, cutoff, prefer = c("above", "below")) {
  prefer <- tryCatch(match.arg(prefer),
                     error = function(e) stop("`prefer` must be \"above\" or \"below\".", call. = FALSE))
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.", call. = FALSE)
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a single non-missing numeric value.", call. = FALSE)
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

#' Prefer individuals with a given value of a covariate
#'
#' @description
#' Scores a candidate set of individuals by the proportion of its members whose
#' `covariate` value equals `level`. Intended for dichotomous covariates, but works for
#' any covariate type supporting `==`.
#'
#' Criteria are used by \code{\link{rkn_select_criteria}}: within each family, where
#' several sets of unrelated individuals remain to choose from, each candidate set is
#' scored and the best one is selected.
#'
#' @param covariate Character. Name of a vertex attribute of the network, such as
#'   \code{"sex"} for \code{V(g)$sex}.
#' @param level The preferred value of `covariate` (matched with `==`).
#'
#' @return An object of class \code{c("rkn_criterion_level", "rkn_criterion")}.
#'
#' @seealso \code{\link{rkn_select_criteria}}, \code{\link{rkn_criterion_cutoff}};
#'   \code{\link{rkn_select_global}} for a target composition of the whole selection.
#'
#' @examples
#' # Prefer women
#' rkn_criterion_level("sex", "F")
#'
#' @export
rkn_criterion_level <- function(covariate, level) {
  if (!is.character(covariate) || length(covariate) != 1L) {
    stop("`covariate` must be a single character string.", call. = FALSE)
  }
  if (length(level) != 1L || is.na(level)) {
    stop("`level` must be a single non-missing value.", call. = FALSE)
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
