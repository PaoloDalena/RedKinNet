#' Plot the full network colored by diagnostic classes (pre-selection classes).
#'
#' This function is graph-first and expects diagnostics to be available from the
#' selection pipeline (i.e., an object produced by rkn_select_random(..., diagnostics=TRUE)).
#' Vertex colors are assigned according to `node_diagnostics$pre_class`.
#'
#' Graphs are assumed to be undirected by design.
#'
#' @param g An igraph object. Vertex identifiers must be stored in `V(g)$name`.
#' @param res A diagnostics-enabled selection output (e.g., from `rkn_select_random()` with
#'   `diagnostics = TRUE`). Must contain `res$node_diagnostics` with columns `id` and `pre_class`.
#' @param layout Character. Graph layout algorithm; one of:
#'   - `"fr"`: Fruchterman-Reingold (default; recommended for multi-component graphs),
#'   - `"gem"`: Graph Embedder Method (may be computationally expensive; can be very slow for
#'     networks with many families/components).
#' @param vertex_size Numeric. Size of vertices in the plot.
#' @param vertex_label Logical. If TRUE, vertex labels (`V(g)$name`) are displayed.
#' @param show_legend Logical. If TRUE, a legend is displayed for the diagnostic classes present.
#' @param main Character. Title of the plot.
#'
#' @return Invisibly returns NULL.
#' @export
rkn_plot_diagnostics <- function(g,
                                 res,
                                 layout = c("fr", "gem"),
                                 vertex_size = 4,
                                 vertex_label = FALSE,
                                 show_legend = TRUE,
                                 main = "Full network (colored by diagnostics)") {

  layout <- match.arg(layout)

  if (!igraph::is_igraph(g)) {
    stop("g must be an igraph object.")
  }
  if (is.null(igraph::V(g)$name)) {
    stop("V(g)$name must be defined (vertex identifiers are required).")
  }
  if (missing(res) || is.null(res$node_diagnostics)) {
    stop("res$node_diagnostics not found. Provide a selection output with diagnostics=TRUE.")
  }

  nd <- res$node_diagnostics
  required_cols <- c("id", "pre_class")
  if (!all(required_cols %in% names(nd))) {
    stop("res$node_diagnostics must contain columns: 'id' and 'pre_class'.")
  }

  vnames <- igraph::V(g)$name
  cls_map <- stats::setNames(rep(NA_character_, length(vnames)), vnames)

  nd$id <- as.character(nd$id)
  nd$pre_class <- as.character(nd$pre_class)

  # If duplicates exist, the last occurrence would overwrite earlier ones.
  # This is acceptable only if duplicates are semantically identical.
  cls_map[intersect(nd$id, vnames)] <- nd$pre_class[match(intersect(nd$id, vnames), nd$id)]

  allowed <- c("kept_apriori", "discarded_apriori", "candidate_simple", "candidate_complex")
  if (any(is.na(cls_map))) {
    stop("Unclassified nodes found: some V(g)$name values have no diagnostic class in res$node_diagnostics.")
  }
  if (any(!cls_map %in% allowed)) {
    bad <- sort(unique(cls_map[!cls_map %in% allowed]))
    stop(paste0(
      "Unexpected diagnostic classes found in res$node_diagnostics$pre_class: ",
      paste(bad, collapse = ", "),
      ". Allowed: ", paste(allowed, collapse = ", "), "."
    ))
  }

  col_map <- c(
    kept_apriori = "#1B9E77",
    discarded_apriori = "#B3B3B3",
    candidate_simple = "#FFD92F",
    candidate_complex = "#E41A1C"
  )
  label_map <- c(
    kept_apriori = "Kept (a priori)",
    discarded_apriori = "Discarded (a priori)",
    candidate_simple = "Candidate (simple)",
    candidate_complex = "Candidate (complex)"
  )

  class_order <- c(
    "kept_apriori",
    "discarded_apriori",
    "candidate_simple",
    "candidate_complex"
  )
  vcol <- unname(col_map[cls_map])

  lay <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    gem = igraph::layout_with_gem(g)
  )

  op <- graphics::par(mar = c(0.5, 0.5, 2, 0.5))
  on.exit(graphics::par(op), add = TRUE)

  igraph::plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) vnames else NA,
    vertex.label.cex = 0.6,
    vertex.label.color = "black",
    edge.arrow.size = 0,   # graphs are undirected by design
    edge.curved = 0.1,
    main = main
  )

  if (isTRUE(show_legend)) {

    present <- class_order[class_order %in% cls_map]

    graphics::legend(
      "topright",
      legend = label_map[present],
      col = unname(col_map[present]),
      pch = 20,
      pt.cex = 1,
      cex = 0.75,
      bty = "o",
      bg = "white",
      box.col = "grey",
      box.lwd = 0.5,
      title = NULL,
      y.intersp = 0.7
    )
  }
}

