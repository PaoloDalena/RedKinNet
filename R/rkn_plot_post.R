#' Plot the full network after selection, highlighting kept vs discarded nodes.
#'
#' This function is graph-first and visualizes the outcome of the selection step.
#' Vertices are colored according to whether they are kept or discarded.
#' Nodes not included in `kept` are treated as discarded.
#'
#' Graphs are assumed to be undirected by design.
#'
#' @param g An igraph object. Vertex identifiers must be stored in `V(g)$name`.
#' @param kept Character vector of vertex identifiers selected by the algorithm.
#' @param layout Character. Graph layout algorithm; one of:
#'   - `"fr"`: Fruchterman-Reingold (default; recommended),
#'   - `"gem"`: Graph Embedder Method (may be slow for networks with many components).
#' @param vertex_size Numeric. Size of vertices in the plot.
#' @param vertex_label Logical. If TRUE, vertex labels (`V(g)$name`) are displayed.
#' @param col_kept Character. Color used for kept vertices (default: "#1B9E77").
#' @param col_discarded Character. Color used for discarded vertices (default: "#B3B3B3").
#' @param main Character. Title of the plot.
#'
#' @return Invisibly returns NULL.
#' @export
rkn_plot_post <- function(g,
                          kept,
                          layout = c("fr", "gem"),
                          vertex_size = 4,
                          vertex_label = FALSE,
                          col_kept = "#1B9E77",
                          col_discarded = "#B3B3B3",
                          main = "Full network post-selection: kept vs discarded") {

  layout <- match.arg(layout)

  if (!igraph::is_igraph(g)) {
    stop("g must be an igraph object.")
  }
  if (is.null(igraph::V(g)$name)) {
    stop("V(g)$name must be defined (vertex identifiers are required).")
  }
  if (missing(kept)) {
    stop("Argument 'kept' must be provided.")
  }

  vnames <- igraph::V(g)$name
  kept <- unique(as.character(kept))
  kept <- kept[kept %in% vnames]

  # ---- Assign colors based on selection outcome -------------------------------
  is_kept <- vnames %in% kept
  vcol <- ifelse(is_kept, col_kept, col_discarded)

  # ---- Compute layout ---------------------------------------------------------
  lay <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    gem = igraph::layout_with_gem(g)
  )

  # ---- Plot -------------------------------------------------------------------
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, 2, 0.5))

  igraph::plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) vnames else NA,
    vertex.label.cex = 0.6,
    edge.arrow.size = 0,   # graphs are undirected by design
    edge.curved = 0.1,
    main = main
  )

  invisible(NULL)
}
