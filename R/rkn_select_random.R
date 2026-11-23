#' Select unrelated individuals at random
#'
#' @description
#' Selects a maximum set of unrelated individuals, a set with no related pair of the
#' largest possible size (see \code{\link{rkn_selection}}), and makes at random the
#' choices that the network leaves open. The network is partitioned into families, its
#' connected components, and each family is resolved on its own: since no relationship
#' crosses two families, the union of their selections has no related pair. By default,
#' each family is drawn uniformly at random among all its maximum sets by an exact
#' algorithm that also handles the intricate families of genetic isolates (see the
#' section Search over maximum sets of \code{\link{rkn_selection}}); alternatively,
#' structural rules decide most individuals.
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
#'   of the returned object: computing the inclusion probabilities makes a selection two
#'   to three times slower in family studies, though hardly slower in population cohorts.
#'   On very large networks, or when many selections are drawn, \code{diagnostics = FALSE}
#'   returns the selection alone.
#' @param uniform Logical. If \code{TRUE} (default), each family is drawn uniformly at
#'   random among all its maximum sets, so that every maximum set of the network has the
#'   same probability, and the selection probabilities are exact: with diagnostics,
#'   \code{node_diagnostics$inclusion_probability} gives the probability that each
#'   individual is selected. Families too intricate for this exact computation are
#'   selected as with \code{uniform = FALSE} and reported in the column \code{uniform}
#'   of \code{family_summary} (see \code{\link{rkn_selection}}). If \code{FALSE},
#'   structural rules decide most individuals and the remaining choices are made at
#'   random; the result is still a maximum set, but the maximum sets are not equally
#'   likely, and the rules never choose some of them: in a chain of six relatives, two
#'   of its four maximum sets. The rules record why each individual was kept or
#'   discarded, and are somewhat faster in cohorts whose families are small.
#'
#' @section Computation time:
#' On one core of a laptop (Apple M2), the uniform selection with diagnostics took about
#' 0.3 seconds for a population cohort of 1,000 samples, 3 seconds for 10,000 and 30
#' seconds for 100,000, and the selection by rules about two thirds of these times. The
#' time depends on the largest and most intricate families more than on the number of
#' samples: a family study of 10,000 samples took 40 seconds (16 without diagnostics),
#' and a family study of 100,000 samples about 12 minutes on one core of a compute node,
#' where a population cohort of the same size took less than 2 minutes.
#'
#' @return An object of class \code{\link{rkn_selection}}, with \code{print},
#'   \code{summary} and \code{plot} methods; \code{$kept} holds the identifiers of the
#'   selected individuals.
#'
#' @seealso \code{\link{rkn_selection}}, \code{\link{summary.rkn_selection}},
#'   \code{\link{plot.rkn_selection}}, \code{\link{rkn_select_criteria}},
#'   \code{\link{rkn_select_global}}.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' # Uniformly among all the maximum sets, with the probability that each individual is
#' # selected
#' sel <- rkn_select_random(g, seed = 1)
#' sel
#' sel$kept
#' sel$node_diagnostics[, c("id", "inclusion_probability")]
#'
#' # The selection alone, without diagnostics
#' rkn_select_random(g, seed = 1, diagnostics = FALSE)$kept
#'
#' # By structural rules, which record why each individual was kept or discarded
#' rules <- rkn_select_random(g, uniform = FALSE, seed = 1)
#' summary(rules)
#'
#' @export
rkn_select_random <- function(g, uniform = TRUE, seed = NULL, diagnostics = TRUE) {

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
