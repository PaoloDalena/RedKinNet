# ---------- Helper: stable palette for categorical family identifiers ----------
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
#' @param layout Character. Graph layout algorithm; one of:
#'   - `"fr"`: Fruchterman-Reingold (default; recommended for multi-component graphs),
#'   - `"gem"`: Graph Embedder Method (can yield aesthetically pleasing layouts but may be
#'     computationally expensive; it can be very slow for networks with many families/components).
#' @param vertex_size Numeric. Size of vertices in the plot.
#' @param vertex_label Logical. If TRUE, vertex labels (`V(g)$name`) are displayed.
#' @param show_legend Logical. If TRUE, a legend mapping colors to family IDs is shown.
#'   For readability, the legend is displayed only if the number of families is ≤ 25.
#' @param main Character. Title of the plot.
#'
#' @return Invisibly returns NULL.
#' @export
rkn_plot_pre <- function(g,
                         layout = c("fr", "gem"),
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

  # ---- Determine family identifiers -------------------------------------------
  # Prefer an existing vertex attribute; otherwise, fall back to connected components.
  fam_ids <- igraph::V(g)$family_id
  if (is.null(fam_ids)) {
    fam_ids <- igraph::components(g)$membership
  }

  fam_levels <- sort(unique(fam_ids))
  pal <- .rkn_palette(length(fam_levels))
  fam_col_map <- stats::setNames(pal, as.character(fam_levels))
  vcol <- unname(fam_col_map[as.character(fam_ids)])

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
    vertex.label = if (isTRUE(vertex_label)) igraph::V(g)$name else NA,
    vertex.label.cex = 0.6,
    edge.arrow.size = 0,   # graphs are undirected by design
    edge.curved = 0.1,
    main = main
  )

  # ---- Legend -----------------------------------------------------------------
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
