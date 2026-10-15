#' Internal helper: plot the network after selection
#'
#' Selected individuals are drawn in the \code{kept} colour of \code{palette}, the
#' others in the \code{discarded} colour. Called by \code{plot.rkn_selection()}, whose
#' documentation describes the arguments.
#'
#' @return The vertex encodings, invisibly.
#'
#' @keywords internal
#' @noRd
.rkn_plot_post <- function(g, kept, palette = "colorblind", layout = "auto", seed = 1,
                           min_family_size = 1, vertex_size = NULL, edge_width = 1,
                           vertex_label = FALSE, show_legend = TRUE, main = NULL) {
  fam <- .rkn_filter_families(g, min_family_size)
  g <- fam$g
  pal <- .rkn_palettes[[palette]]
  is_kept <- igraph::V(g)$name %in% kept
  n <- length(is_kept)

  if (is.null(main)) {
    main <- sprintf("Selected %s of %s individuals%s (%.1f%%)",
                    format(sum(is_kept), big.mark = ","), format(n, big.mark = ","),
                    if (min_family_size > 1) sprintf(" in families of %d+", min_family_size) else "",
                    100 * sum(is_kept) / n)
  }
  legend <- list(labels = c("selected", "discarded"),
                 fill = unname(pal[c("kept", "discarded")]), ncol = 2)

  .rkn_draw_network(g, fill = unname(ifelse(is_kept, pal["kept"], pal["discarded"])),
                    frame_col = "grey25", frame_width = 1,
                    layout = layout, seed = seed, family_id = fam$family_id,
                    vertex_size = vertex_size,
                    edge_width = edge_width, vertex_label = vertex_label,
                    legend = if (isTRUE(show_legend)) legend, main = main)
}
