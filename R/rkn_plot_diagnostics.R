#' Internal helper: plot the network by type of decision
#'
#' Vertices are coloured by \code{res$node_diagnostics$pre_class}, the way each
#' individual was decided. With \code{show_outcome = TRUE}, selected individuals are
#' filled and discarded individuals are hollow, with a frame in the colour of their
#' class. Called by \code{plot.rkn_selection()}, whose documentation describes the
#' arguments.
#'
#' @return The vertex encodings, invisibly.
#'
#' @keywords internal
#' @noRd
.rkn_plot_diagnostics <- function(g, res, palette = "colorblind", show_outcome = TRUE,
                                  layout = "auto", seed = 1, vertex_size = NULL,
                                  edge_width = 1, vertex_label = FALSE,
                                  show_legend = TRUE, main = NULL) {
  pal <- .rkn_palettes[[palette]]
  nd <- res$node_diagnostics
  row <- match(igraph::V(g)$name, nd$id)
  decision <- nd$pre_class[row]
  is_kept <- nd$post_class[row] == "kept"
  n <- length(is_kept)

  fill <- unname(pal[decision])
  frame_col <- rep("grey25", n)
  frame_width <- rep(1, n)
  if (isTRUE(show_outcome)) {
    frame_col[!is_kept] <- fill[!is_kept]
    fill[!is_kept] <- "white"
    frame_width[!is_kept] <- 2
  }

  if (is.null(main)) {
    main <- sprintf("How each individual was decided (%s of %s selected)",
                    format(sum(is_kept), big.mark = ","), format(n, big.mark = ","))
  }
  present <- names(.rkn_decision_labels)[names(.rkn_decision_labels) %in% decision]
  legend <- list(labels = unname(.rkn_decision_labels[present]), fill = unname(pal[present]),
                 ncol = 2,
                 title = if (isTRUE(show_outcome)) "Filled: selected. Hollow: discarded.")

  .rkn_draw_network(g, fill = fill, frame_col = frame_col, frame_width = frame_width,
                    layout = layout, seed = seed, vertex_size = vertex_size,
                    edge_width = edge_width, vertex_label = vertex_label,
                    legend = if (isTRUE(show_legend)) legend, main = main)
}
