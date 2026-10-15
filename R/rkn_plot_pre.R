#' Plot a relatedness network
#'
#' @description
#' Draws a relatedness network in which each family (connected component) is laid out
#' on its own.
#' \itemize{
#'   \item \code{rkn_plot_pre()} draws any network, for instance before a selection,
#'     with the vertices coloured by the size of their family or by family.
#'   \item \code{plot()} draws a selection returned by \code{\link{rkn_select_random}}
#'     or \code{\link{rkn_select_criteria}} (an \code{\link{rkn_selection}} object): the
#'     selected and discarded individuals, how each individual was decided, or the
#'     families before the selection, as \code{rkn_plot_pre()} does.
#' }
#' All plots of the same network put each individual in the same position, so they can
#' be compared or combined in panels, for instance with \code{par(mfrow = c(1, 3))}.
#'
#' @param x,g The network to draw: for \code{plot()}, an \code{\link{rkn_selection}}
#'   object; for \code{rkn_plot_pre()}, an igraph object with the vertex identifiers
#'   in \code{V(g)$name}.
#' @param type What the colours show (\code{plot()} only):
#'   \itemize{
#'     \item \code{"post"} (default): selected and discarded individuals;
#'     \item \code{"diagnostics"}: how each individual was decided (see Details);
#'       requires a selection computed with \code{diagnostics = TRUE};
#'     \item \code{"pre"}: the families before the selection, as drawn by
#'       \code{rkn_plot_pre()}.
#'   }
#' @param palette Character. Colours of the \code{"post"} and \code{"diagnostics"}
#'   plots (\code{plot()} only). \code{"colorblind"} (default) uses the Okabe-Ito
#'   colours, which remain distinguishable with the most common forms of colour
#'   blindness; \code{"classic"} uses the earlier colours of the package, which pair red
#'   with green.
#' @param show_outcome Logical. In the \code{"diagnostics"} plot (\code{plot()} only),
#'   if \code{TRUE} (default), selected individuals are drawn as filled vertices and
#'   discarded individuals as hollow vertices with a frame in the colour of their
#'   class, so the plot shows both how and with which outcome each individual was
#'   decided.
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
#' @param colour_by Character. \code{"size"} (default) colours each vertex by the size
#'   of its family, in the classes 1, 2, 3, 4, 5 and 6 or more; \code{"family"} gives
#'   each family its own colour, which is useful only for a handful of families and is
#'   drawn without a legend.
#' @param vertex_size Numeric. Size of the vertices. By default (\code{NULL}) it
#'   decreases with the number of individuals. Enlarging the graphics device does not
#'   separate overlapping vertices, because vertex sizes are relative to the plotting
#'   region; reduce \code{vertex_size} instead.
#' @param edge_width Numeric. Width of the edges.
#' @param vertex_label Logical. If \code{TRUE}, the vertex identifiers are written next
#'   to the vertices; this is readable only for small networks.
#' @param show_legend Logical. If \code{TRUE} (default), a legend is drawn below the
#'   network.
#' @param main Character. Title of the plot. By default (\code{NULL}) an informative
#'   title is generated, such as the number of individuals and families or of selected
#'   individuals; replace it with \code{main = "My title"} or remove it with
#'   \code{main = ""}.
#'
#' @return \code{rkn_plot_pre()} returns \code{NULL} and \code{plot()} returns \code{x},
#'   both invisibly.
#'
#' @seealso \code{\link{rkn_selection}}
#'
#' @order 2
#' @export
rkn_plot_pre <- function(g,
                         layout = c("auto", "components", "grid"),
                         seed = 1,
                         colour_by = c("size", "family"),
                         vertex_size = NULL,
                         edge_width = 1,
                         vertex_label = FALSE,
                         show_legend = TRUE,
                         main = NULL) {

  layout <- match.arg(layout)
  colour_by <- match.arg(colour_by)

  if (!igraph::is_igraph(g)) {
    stop("g must be an igraph object.")
  }
  if (is.null(igraph::V(g)$name)) {
    stop("V(g)$name must be defined (vertex identifiers are required).")
  }

  comp <- igraph::components(g)
  fam_size <- comp$csize[comp$membership]

  if (identical(colour_by, "size")) {
    size_class <- ifelse(fam_size >= 6, "6+", as.character(fam_size))
    fill <- unname(.rkn_size_class_colours[size_class])
    present <- names(.rkn_size_class_colours)[names(.rkn_size_class_colours) %in% size_class]
    legend <- list(labels = present, fill = unname(.rkn_size_class_colours[present]),
                   ncol = length(present), title = "Family size")
  } else {
    hues <- seq(15, 375, length.out = comp$no + 1)[seq_len(comp$no)]
    fill <- grDevices::hcl(h = hues, c = 70, l = 60)[comp$membership]
    legend <- NULL
  }

  if (is.null(main)) {
    main <- sprintf("%s individuals in %s families",
                    format(igraph::vcount(g), big.mark = ","), format(comp$no, big.mark = ","))
  }

  .rkn_draw_network(g, fill = fill, frame_col = "grey25", frame_width = 1,
                    layout = layout, seed = seed, vertex_size = vertex_size,
                    edge_width = edge_width, vertex_label = vertex_label,
                    legend = if (isTRUE(show_legend)) legend, main = main)
  invisible(NULL)
}
