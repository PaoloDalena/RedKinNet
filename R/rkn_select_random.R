#' Random selection over an entire network (connected-component-wise)
#'
#' @description
#' The network is partitioned into families, defined as the connected components of
#' \code{g}, and each family is resolved independently. Since no relationship crosses
#' two components, the union of the per-family selections is a set of mutually
#' unrelated individuals. Where several maximum sets of unrelated individuals are
#' possible, the choice among them is made at random. Structural rules decide most
#' individuals; those left undecided in large families are resolved by an exact
#' algorithm that also handles the intricate families of genetic isolates (see the
#' section Search over maximum sets of \code{\link{rkn_selection}}).
#'
#' @param g An igraph object. Vertex names must be stored in V(g)$name. Vertex
#'   attributes cannot be named \code{degree}, \code{max_path_length} or \code{.name}:
#'   these names are reserved for the structural measures computed during the selection.
#'   The network is usually built with \code{\link{rkn_graph_from_kinship}}, which also
#'   checks the samples; for a network built in another way, see
#'   \code{\link{rkn_check_samples}}.
#' @param seed Optional integer. If provided, the selection is reproducible: a seed is
#'   drawn from it for each family, so that the selection within a family does not
#'   depend on the random numbers used by the others, and different values give
#'   independent selections. If \code{NULL} (default), the current random number stream
#'   is used.
#' @param diagnostics Logical. If \code{TRUE} (default), the returned object also records
#'   how each individual and each family was resolved (\code{node_diagnostics} and
#'   \code{family_summary}), which \code{print()}, \code{summary()} and
#'   \code{plot(type = "diagnostics")} use. Diagnostics increase the run time and the size
#'   of the returned object; on very large networks, \code{diagnostics = FALSE} returns
#'   the selection alone.
#' @param uniform Logical. If \code{FALSE} (default), structural rules decide most
#'   individuals and the remaining choices are made at random; the result is a maximum
#'   set of unrelated individuals, but not every maximum set has the same probability of
#'   being chosen. If \code{TRUE}, each family is drawn uniformly at random among all
#'   its maximum sets, which makes the selection probabilities exact: with diagnostics,
#'   \code{node_diagnostics$inclusion_probability} gives the probability that each
#'   individual is selected. Families too intricate for this exact computation are
#'   selected as with \code{uniform = FALSE} and reported in the column \code{uniform}
#'   of \code{family_summary} (see \code{\link{rkn_selection}}).
#'
#' @return An object of class \code{\link{rkn_selection}}, with \code{print},
#'   \code{summary} and \code{plot} methods; \code{$kept} holds the identifiers of the
#'   selected individuals.
#'
#' @seealso \code{\link{rkn_selection}}, \code{\link{summary.rkn_selection}},
#'   \code{\link{plot.rkn_selection}}, \code{\link{rkn_select_criteria}}.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' sel <- rkn_select_random(g, seed = 1)
#' sel
#' sel$kept
#' summary(sel)
#'
#' # The selection alone, without diagnostics
#' rkn_select_random(g, seed = 1, diagnostics = FALSE)$kept
#'
#' @export
rkn_select_random <- function(g, uniform = FALSE, seed = NULL, diagnostics = TRUE) {

  .rkn_check_network(g)
  .rkn_check_seed(seed)
  .rkn_check_flag(diagnostics, "diagnostics")
  .rkn_check_flag(uniform, "uniform")

  select_upto5 <- .rkn_select_random_upto5
  select_6plus <- .rkn_select_random_6plus
  if (uniform) {
    select_upto5 <- function(fam_df, g_sub, seed, return_diagnostics) {
      .rkn_select_uniform_family(fam_df, g_sub, seed, return_diagnostics)
    }
    select_6plus <- select_upto5
  }

  .rkn_select_network(
    g,
    select_upto5 = select_upto5,
    select_6plus = select_6plus,
    seed = seed,
    diagnostics = diagnostics,
    method = "random",
    criteria = NULL,
    call = match.call(),
    uniform = uniform
  )
}
