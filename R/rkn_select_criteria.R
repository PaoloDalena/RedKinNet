#' Select unrelated individuals by criteria
#'
#' @description
#' Selects a maximum set of unrelated individuals, as \code{\link{rkn_select_random}}
#' does, but makes the choices that the network leaves open by criteria on the
#' covariates instead of at random. Each family is resolved on its own: families of up to
#' 5 members through the structural rules in \code{\link{rules_db}}, larger families
#' through structural reductions followed by a search of what remains; the returned
#' object is as for \code{rkn_select_random()}. How the criteria are evaluated in large
#' families, exactly or on a sample of maximum sets, is described in the section Search
#' over maximum sets of \code{\link{rkn_selection}}. With \code{optimal = TRUE}, the
#' structural rules are skipped and the criteria choose among all the maximum sets of
#' each family.
#'
#' Covariates referenced by the criteria are read from the vertex attributes of
#' \code{g}. The criteria express preferences for individuals, such as women or people
#' younger than 50; for a target composition of the whole selection, such as half women
#' or a mean age, see \code{\link{rkn_select_global}}.
#'
#' @inheritParams rkn_select_random
#' @param criterion_primary A criterion, created by \code{\link{rkn_criterion_level}} or
#'   \code{\link{rkn_criterion_cutoff}}, used to score candidate resolutions within each
#'   family.
#' @param criterion_secondary Optional \code{rkn_criterion} object used to break ties in
#'   \code{criterion_primary}.
#' @param seed Optional integer, as in \code{\link{rkn_select_random}}: a seed is drawn
#'   from it for each family. Random numbers are used only to break exact ties in the
#'   criterion scores.
#' @param optimal Logical. If \code{FALSE} (default), structural rules decide part of
#'   each family before the criteria are applied, and some of them fix individuals
#'   regardless of the criteria (for instance, they keep the individual at the end of a
#'   chain of relatives); the result is the best maximum set among those compatible with
#'   the rules. If \code{TRUE}, the criteria are evaluated over all the maximum sets of
#'   each family, exactly, or as described in the section Search over maximum sets of
#'   \code{\link{rkn_selection}} when the covariate has missing values. Families too
#'   intricate for this
#'   computation are selected as with \code{optimal = FALSE} and reported in the column
#'   \code{optimal} of \code{family_summary}.
#'
#' @section Computation time:
#' On one core of a laptop (Apple M2), the selection took about as long as
#' \code{rkn_select_random(uniform = FALSE)}: 1.6 seconds for a population cohort of
#' 10,000 samples, 17 seconds for 100,000 and 32 seconds for a family study of 10,000;
#' with \code{optimal = TRUE}, 5.5 seconds, 1 minute and 67 seconds.
#'
#' @return An object of class \code{\link{rkn_selection}}, as for
#'   \code{\link{rkn_select_random}}; its \code{criteria} element holds the criteria.
#'
#' @seealso \code{\link{rkn_selection}}, \code{\link{summary.rkn_selection}},
#'   \code{\link{plot.rkn_selection}}, \code{\link{rkn_select_random}},
#'   \code{\link{rkn_select_global}}, \code{\link{rkn_criterion_level}},
#'   \code{\link{rkn_criterion_cutoff}}, \code{\link{rules_db}}.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' # Where the structure of a family leaves a choice, prefer women ...
#' rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1)$kept
#'
#' # ... evaluating the criterion over all the maximum sets of each family ...
#' rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1, optimal = TRUE)$kept
#'
#' # ... or prefer individuals younger than 40
#' rkn_select_criteria(g, rkn_criterion_cutoff("age", 40, prefer = "below"), seed = 1)$kept
#'
#' # A secondary criterion breaks the ties left by the primary one
#' sel <- rkn_select_criteria(g, rkn_criterion_level("sex", "F"),
#'                            rkn_criterion_cutoff("age", 40, prefer = "below"), seed = 1)
#' sel
#'
#' @export
rkn_select_criteria <- function(g,
                                criterion_primary,
                                criterion_secondary = NULL,
                                optimal = FALSE,
                                seed = NULL,
                                diagnostics = TRUE) {

  .rkn_check_network(g)
  if (missing(criterion_primary) || !inherits(criterion_primary, "rkn_criterion")) {
    stop("`criterion_primary` must be a criterion, created by rkn_criterion_level() or rkn_criterion_cutoff().",
         call. = FALSE)
  }
  if (!is.null(criterion_secondary) && !inherits(criterion_secondary, "rkn_criterion")) {
    stop("`criterion_secondary` must be NULL or a criterion, created by rkn_criterion_level() or ",
         "rkn_criterion_cutoff().", call. = FALSE)
  }
  # The covariates are checked here, since they are read only in the families that
  # leave a choice.
  for (criterion in list(criterion_primary, criterion_secondary)) {
    if (is.null(criterion)) next
    covariate <- criterion$params$covariate
    if (!covariate %in% igraph::vertex_attr_names(g)) {
      stop(sprintf("The network has no vertex attribute \"%s\".", covariate), call. = FALSE)
    }
    if (identical(criterion$type, "cutoff") && !is.numeric(igraph::vertex_attr(g, covariate))) {
      stop(sprintf("A cutoff criterion needs a numeric covariate, and \"%s\" is not numeric.", covariate),
           call. = FALSE)
    }
  }
  .rkn_check_seed(seed)
  .rkn_check_flag(diagnostics, "diagnostics")
  .rkn_check_flag(optimal, "optimal")

  select_upto5 <- function(fam_df, g_sub, seed, return_diagnostics) {
    .rkn_select_criteria_upto5(
      fam_df = fam_df,
      criterion_primary = criterion_primary,
      criterion_secondary = criterion_secondary,
      seed = seed,
      return_diagnostics = return_diagnostics,
      g_sub = g_sub
    )
  }
  select_6plus <- function(fam_df, g_sub, seed, return_diagnostics) {
    .rkn_select_criteria_6plus(
      fam_df = fam_df,
      g_sub = g_sub,
      criterion_primary = criterion_primary,
      criterion_secondary = criterion_secondary,
      seed = seed,
      return_diagnostics = return_diagnostics
    )
  }
  if (optimal) {
    select_upto5 <- function(fam_df, g_sub, seed, return_diagnostics) {
      .rkn_select_optimal_family(fam_df, g_sub, criterion_primary, criterion_secondary, seed, return_diagnostics)
    }
    select_6plus <- select_upto5
  }

  .rkn_select_network(
    g,
    select_upto5 = select_upto5,
    select_6plus = select_6plus,
    seed = seed,
    diagnostics = diagnostics,
    method = "criteria",
    criteria = list(primary = criterion_primary, secondary = criterion_secondary),
    call = match.call(),
    optimal = optimal
  )
}
