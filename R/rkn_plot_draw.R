.rkn_palettes <- list(
  colorblind = c(kept = "#009E73", discarded = "#999999",
                 kept_apriori = "#009E73", discarded_apriori = "#999999",
                 candidate_simple = "#E69F00", candidate_complex = "#0072B2"),
  classic = c(kept = "#1B9E77", discarded = "#B3B3B3",
              kept_apriori = "#1B9E77", discarded_apriori = "#B3B3B3",
              candidate_simple = "#FFD92F", candidate_complex = "#E41A1C")
)

.rkn_size_class_colours <- c("1" = "#C6DBEF", "2" = "#9ECAE1", "3" = "#6BAED6",
                             "4" = "#4292C6", "5" = "#2171B5", "6+" = "#084594")

.rkn_decision_labels <- c(kept_apriori = "kept by a structural rule",
                          discarded_apriori = "discarded by a structural rule",
                          candidate_simple = "decided among single candidates",
                          candidate_complex = "decided among alternative sets")

#' Internal helper: keep the families with at least a given number of members
#'
#' Returns the subgraph of the families of \code{g} with at least
#' \code{min_family_size} members and, for each of its vertices, the identifier of its
#' family in \code{g}, which the layout uses to keep the shape of each family.
#'
#' @keywords internal
#' @noRd
.rkn_filter_families <- function(g, min_family_size) {
  comp <- igraph::components(g)
  keep <- which(comp$csize[comp$membership] >= min_family_size)
  if (length(keep) == 0) stop("No family has at least ", min_family_size, " members.", call. = FALSE)
  list(g = igraph::induced_subgraph(g, keep), family_id = comp$membership[keep])
}

#' Internal helper: draw a network of families
#'
#' Lays out \code{g} with \code{.rkn_layout()}, passing on \code{family_id}, and draws it
#' with straight edges. When
#' \code{vertex_size} is \code{NULL}, it decreases with the number of vertices, and
#' frames are thinned on vertices smaller than size 4. The
#' legend, a list with elements \code{labels}, \code{fill}, \code{ncol} and optionally
#' \code{title}, is drawn below the network, in a bottom margin sized to its number of
#' rows, so that it never covers a vertex. Only the margins are changed, and they are
#' restored on exit.
#'
#' @return The vertex encodings (\code{fill}, \code{frame_col}, \code{frame_width}),
#'   invisibly.
#'
#' @keywords internal
#' @noRd
.rkn_draw_network <- function(g, fill, frame_col, frame_width, layout, seed, family_id,
                              vertex_size, edge_width, vertex_label, legend, main) {
  xy <- .rkn_layout(g, layout, seed, family_id)
  if (is.null(vertex_size)) vertex_size <- min(6, max(0.5, 55 / sqrt(igraph::vcount(g))))
  # On small vertices a full-width frame would hide the fill colour.
  frame_width <- frame_width * min(1, vertex_size / 4)

  legend_lines <- 0
  if (!is.null(legend)) {
    legend_lines <- 0.9 * (ceiling(length(legend$labels) / legend$ncol) + !is.null(legend$title)) + 0.3
  }
  op <- graphics::par(mar = c(0.5 + legend_lines, 0.5, if (nzchar(main)) 2 else 0.5, 0.5))
  on.exit(graphics::par(op), add = TRUE)

  igraph::plot.igraph(
    g,
    layout = xy,
    rescale = FALSE,
    xlim = c(-1.05, 1.05),
    ylim = c(-1.05, 1.05),
    vertex.color = fill,
    vertex.frame.color = frame_col,
    vertex.frame.width = frame_width,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) igraph::V(g)$name else NA,
    vertex.label.cex = 0.55,
    vertex.label.dist = 0.9,
    vertex.label.degree = -pi / 4,
    vertex.label.color = "grey15",
    vertex.label.family = "sans",
    edge.color = "grey45",
    edge.width = edge_width,
    edge.curved = 0,
    edge.arrow.size = 0,
    main = main
  )

  if (!is.null(legend)) {
    usr <- graphics::par("usr")
    graphics::legend(
      x = mean(usr[1:2]), y = usr[3], xjust = 0.5, yjust = 1, xpd = NA,
      legend = legend$labels, title = legend$title, ncol = legend$ncol,
      pch = 21, pt.bg = legend$fill, col = "grey25", pt.cex = 1.3,
      cex = 0.8, bty = "n"
    )
  }
  invisible(list(fill = fill, frame_col = frame_col, frame_width = frame_width))
}
