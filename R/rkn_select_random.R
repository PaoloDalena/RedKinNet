#' Random selection over an entire network (connected-component-wise)
#'
#' @description
#' The network is partitioned into families, defined as the connected components of
#' \code{g}, and each family is resolved independently. Since no relationship crosses
#' two components, the union of the per-family selections is a set of mutually
#' unrelated individuals.
#'
#' @param g An igraph object. Vertex names must be stored in V(g)$name.
#' @param seed Optional integer base seed. If provided, family-specific seeds are derived as seed + family_id.
#' @param diagnostics Logical. If TRUE, returns node-level diagnostics and a standardized family summary.
#' @param return_results_by_family Logical. If TRUE (and diagnostics=TRUE), also returns the per-family raw outputs
#'   in a list (potentially large object).
#' @param verbose Logical. If TRUE, prints a formatted summary to the console.
#' @param plot Character. One of:
#'   \itemize{
#'     \item \code{"no"}: do not produce any plot (default);
#'     \item \code{"pre"}: plot the full network pre-selection, colored by family (calls \code{rkn_plot_pre()});
#'     \item \code{"post"}: plot the full network post-selection, colored by kept vs discarded;
#'     \item \code{"diagnostics"}: plot the full network colored by diagnostic classes.
#'   }
#'   If \code{plot="diagnostics"} is requested with \code{diagnostics=FALSE}, a warning is issued and the function
#'   falls back to \code{plot="post"} (kept vs discarded).
#' @param plot_args Optional named list of arguments passed to the plotting function; see
#'   \code{\link{rkn_plot_pre}} and \code{\link{plot.rkn_selection}}.
#'
#' @return If diagnostics=FALSE: a character vector of selected node IDs (kept).
#'         If diagnostics=TRUE: a list with elements:
#'           kept, discarded, node_diagnostics, family_summary,
#'           and optionally results_by_family (if return_results_by_family=TRUE).
#'
#' @seealso \code{\link{rkn_plot_pre}}, \code{\link{plot.rkn_selection}}.
#'
#' @export
rkn_select_random <- function(g,
                              seed = NULL,
                              diagnostics = FALSE,
                              return_results_by_family = FALSE,
                              verbose = TRUE,
                              plot = c("no", "pre", "post", "diagnostics"),
                              plot_args = list()) {

  plot <- match.arg(plot)

  if (!igraph::is_igraph(g)) stop("g must be an igraph object.")
  if (igraph::vcount(g) == 0) stop("g must contain at least one vertex.")
  if (is.null(igraph::V(g)$name)) stop("V(g)$name must be defined (vertex identifiers are required).")
  if (!is.list(plot_args)) stop("plot_args must be a named list.")

  .rkn_select_network(
    g,
    select_upto5 = rkn_select_random_upto5,
    select_6plus = rkn_select_random_6plus,
    seed = seed,
    diagnostics = diagnostics,
    return_results_by_family = return_results_by_family,
    verbose = verbose,
    plot = plot,
    plot_args = plot_args,
    caller = "rkn_select_random"
  )
}
