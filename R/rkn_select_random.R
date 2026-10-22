#' Random selection over an entire network (connected-component-wise)
#'
#' @description
#' The network is partitioned into families, defined as the connected components of
#' \code{g}, and each family is resolved independently. Since no relationship crosses
#' two components, the union of the per-family selections is a set of mutually
#' unrelated individuals. Where several maximum sets of unrelated individuals are
#' possible, the choice among them is made at random.
#'
#' @param g An igraph object. Vertex names must be stored in V(g)$name. Vertex
#'   attributes cannot be named \code{degree}, \code{max_path_length} or \code{.name}:
#'   these names are reserved for the structural measures computed during the selection.
#'   The network is usually built with \code{\link{rkn_graph_from_kinship}}, which also
#'   checks the samples; for a network built in another way, see
#'   \code{\link{rkn_check_samples}}.
#' @param seed Optional integer base seed. If provided, family-specific seeds are derived as seed + family_id.
#' @param diagnostics Logical. If \code{TRUE} (default), the returned object also records
#'   how each individual and each family was resolved (\code{node_diagnostics} and
#'   \code{family_summary}), which \code{print()}, \code{summary()} and
#'   \code{plot(type = "diagnostics")} use. Diagnostics increase the run time and the size
#'   of the returned object; on very large networks, \code{diagnostics = FALSE} returns
#'   the selection alone.
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
rkn_select_random <- function(g, seed = NULL, diagnostics = TRUE) {

  if (!igraph::is_igraph(g)) stop("g must be an igraph object.")
  if (igraph::vcount(g) == 0) stop("g must contain at least one vertex.")
  if (is.null(igraph::V(g)$name)) stop("V(g)$name must be defined (vertex identifiers are required).")
  reserved <- intersect(igraph::vertex_attr_names(g), .rkn_reserved_attributes)
  if (length(reserved) > 0) stop(.rkn_reserved_message(reserved))

  .rkn_select_network(
    g,
    select_upto5 = rkn_select_random_upto5,
    select_6plus = rkn_select_random_6plus,
    seed = seed,
    diagnostics = diagnostics,
    method = "random",
    criteria = NULL,
    call = match.call()
  )
}
