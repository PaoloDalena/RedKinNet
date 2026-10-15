.rkn_palette <- function(n) {
  grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[1:n],
    c = 70,
    l = 55
  )
}

#' Plot the full network before selection, coloring vertices by family.
#'
#' This function is graph-first and does not depend on the selection step.
#' Vertices are colored according to their family membership. If the vertex
#' attribute `family_id` is available, it is used directly; otherwise,
#' connected components are computed internally and used as family identifiers.
#'
#' Graphs are assumed to be undirected by design.
#'
#' @param g An igraph object. Vertex identifiers must be stored in `V(g)$name`.
#' @param layout Character. How the network is laid out. Each family (connected
#'   component) is always drawn on its own; the options differ in how the families are
#'   arranged:
#'   \itemize{
#'     \item \code{"auto"} (default): \code{"components"} for networks of up to 500
#'       individuals, \code{"grid"} for larger networks;
#'     \item \code{"components"}: the families are packed together around the centre;
#'     \item \code{"grid"}: each family gets its own cell, ordered by family size
#'       (largest first), which keeps large networks readable.
#'   }
#' @param seed Integer. Seed of the random layout, applied locally: the random number
#'   generator of the session is left untouched. The same seed always gives the same
#'   coordinates, so plots of the same network can be compared; another value gives a
#'   different arrangement.
#' @param vertex_size Numeric. Size of vertices in the plot.
#' @param vertex_label Logical. If TRUE, vertex labels (`V(g)$name`) are displayed.
#' @param show_legend Logical. If TRUE, a legend mapping colors to family IDs is shown.
#'   For readability, the legend is displayed only if there are at most 25 families.
#' @param main Character. Title of the plot.
#'
#' @return Invisibly returns NULL.
#' @export
rkn_plot_pre <- function(g,
                         layout = c("auto", "components", "grid"),
                         seed = 1,
                         vertex_size = 4,
                         vertex_label = FALSE,
                         show_legend = FALSE,
                         main = "Full network pre-selection: families") {

  layout <- match.arg(layout)

  if (!igraph::is_igraph(g)) {
    stop("g must be an igraph object.")
  }
  if (is.null(igraph::V(g)$name)) {
    stop("V(g)$name must be defined (vertex identifiers are required).")
  }

  # Prefer an existing vertex attribute; otherwise, fall back to connected components.
  fam_ids <- igraph::V(g)$family_id
  if (is.null(fam_ids)) {
    fam_ids <- igraph::components(g)$membership
  }

  fam_levels <- sort(unique(fam_ids))
  pal <- .rkn_palette(length(fam_levels))
  fam_col_map <- stats::setNames(pal, as.character(fam_levels))
  vcol <- unname(fam_col_map[as.character(fam_ids)])

  lay <- .rkn_layout(g, layout, seed)

  op <- graphics::par(mar = c(0.5, 0.5, 2, 0.5))
  on.exit(graphics::par(op), add = TRUE)

  igraph::plot.igraph(
    g,
    layout = lay,
    rescale = FALSE,
    xlim = c(-1.05, 1.05),
    ylim = c(-1.05, 1.05),
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) igraph::V(g)$name else NA,
    vertex.label.cex = 0.6,
    edge.arrow.size = 0,   # graphs are undirected by design
    edge.curved = 0.1,
    main = main
  )

  if (isTRUE(show_legend)) {
    if (length(fam_levels) <= 25) {
      graphics::legend(
        "topleft",
        legend = paste0("F", fam_levels),
        col = pal,
        pch = 19,
        cex = 0.7,
        bty = "n",
        title = "family_id"
      )
    } else {
      graphics::mtext(
        side = 3,
        line = 0.2,
        adj = 1,
        cex = 0.8,
        paste0(
          "Legend hidden (", length(fam_levels),
          " families). Set show_legend=FALSE or reduce the number of families."
        )
      )
    }
  }

  invisible(NULL)
}
