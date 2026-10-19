#' Criteria-based selection over an entire network (connected-component-wise)
#'
#' @description
#' Criteria-based counterpart of \code{\link{rkn_select_random}}. The network is
#' partitioned into families, defined as the connected components of \code{g}, and
#' each family is resolved independently: families with at most 5 nodes through the
#' structural rules in \code{\link{rules_db}}, larger families through structural
#' reductions followed by the resolution of the remaining ambiguity. Ambiguous
#' selections are resolved by criterion score rather than at random; the partition into
#' families and the returned object are as in \code{rkn_select_random()}.
#'
#' Covariates referenced by the criteria are read from the vertex attributes of
#' \code{g}.
#'
#' @inheritParams rkn_select_random
#' @param criterion_primary An \code{rkn_criterion} object (see
#'   \code{\link{rkn_criterion_target_mean}} and siblings) used to score candidate
#'   resolutions within each family.
#' @param criterion_secondary Optional \code{rkn_criterion} object used to break ties in
#'   \code{criterion_primary}.
#' @param seed Optional integer base seed. If provided, family-specific seeds are derived
#'   as seed + family_id and are used only to break exact ties in the criterion score(s).
#'
#' @return An object of class \code{\link{rkn_selection}}, as for
#'   \code{\link{rkn_select_random}}; its \code{criteria} element holds the criteria.
#'
#' @seealso \code{\link{rkn_selection}}, \code{\link{summary.rkn_selection}},
#'   \code{\link{plot.rkn_selection}}, \code{\link{rkn_select_random}},
#'   \code{\link{rkn_criterion_target_mean}}, \code{\link{rules_db}}.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' # Where the structure of a family leaves a choice, prefer women ...
#' rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1)$kept
#'
#' # ... or a mean age close to 40
#' rkn_select_criteria(g, rkn_criterion_target_mean("age", target = 40, tolerance = 20),
#'                     seed = 1)$kept
#'
#' # A secondary criterion breaks the ties left by the primary one
#' sel <- rkn_select_criteria(g, rkn_criterion_level("sex", "F"),
#'                            rkn_criterion_target_mean("age", target = 40, tolerance = 20),
#'                            seed = 1)
#' sel
#'
#' @export
rkn_select_criteria <- function(g,
                                criterion_primary,
                                criterion_secondary = NULL,
                                seed = NULL,
                                diagnostics = TRUE) {

  if (!igraph::is_igraph(g)) stop("g must be an igraph object.")
  if (igraph::vcount(g) == 0) stop("g must contain at least one vertex.")
  if (is.null(igraph::V(g)$name)) stop("V(g)$name must be defined (vertex identifiers are required).")
  if (missing(criterion_primary) || !inherits(criterion_primary, "rkn_criterion")) {
    stop("`criterion_primary` must be an rkn_criterion object (see rkn_criterion_*()).")
  }
  if (!is.null(criterion_secondary) && !inherits(criterion_secondary, "rkn_criterion")) {
    stop("`criterion_secondary` must be an rkn_criterion object or NULL.")
  }

  .rkn_select_network(
    g,
    select_upto5 = function(fam_df, g_sub, seed, return_diagnostics) {
      rkn_select_criteria_upto5(
        fam_df = fam_df,
        criterion_primary = criterion_primary,
        criterion_secondary = criterion_secondary,
        seed = seed,
        return_diagnostics = return_diagnostics,
        g_sub = g_sub
      )
    },
    select_6plus = function(fam_df, g_sub, seed, return_diagnostics) {
      rkn_select_criteria_6plus(
        fam_df = fam_df,
        g_sub = g_sub,
        criterion_primary = criterion_primary,
        criterion_secondary = criterion_secondary,
        seed = seed,
        return_diagnostics = return_diagnostics
      )
    },
    seed = seed,
    diagnostics = diagnostics,
    method = "criteria",
    criteria = list(primary = criterion_primary, secondary = criterion_secondary),
    call = match.call()
  )
}
