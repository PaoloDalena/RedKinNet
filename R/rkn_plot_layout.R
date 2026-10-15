#' Internal helper: coordinates of a network of separate families
#'
#' Each family (connected component) is laid out on its own with the
#' Fruchterman-Reingold algorithm, using the seed \code{seed + k} for the family with
#' identifier k. Identifiers are the component numbers of \code{g} unless
#' \code{family_id} gives them, one per vertex: a subgraph of whole families then keeps
#' the identifiers, hence the shapes, of the network it was taken from.
#' The families are then arranged on the page: packed together by diffusion-limited
#' aggregation (\code{"components"}, as in \code{igraph::layout_components()}) or given
#' one cell each, by decreasing size (\code{"grid"}); \code{"auto"} chooses the former
#' up to 500 vertices and the latter above. The caller's random number stream is left
#' untouched. Coordinates are centred and scaled between -1 and 1 by a common factor, which
#' preserves the aspect ratio.
#'
#' @keywords internal
#' @noRd
.rkn_layout <- function(g, layout, seed, family_id = NULL) {
  n <- igraph::vcount(g)
  if (identical(layout, "auto")) layout <- if (n <= 500) "components" else "grid"
  withr::local_preserve_seed()

  if (is.null(family_id)) family_id <- igraph::components(g)$membership
  fams <- split(seq_len(n), family_id)
  ids <- as.integer(names(fams))
  subs <- vector("list", length(fams))
  local <- vector("list", length(fams))
  for (k in seq_along(fams)) {
    subs[[k]] <- igraph::induced_subgraph(g, fams[[k]])
    if (length(fams[[k]]) == 1) {
      local[[k]] <- matrix(0, 1, 2)
    } else {
      set.seed(seed + ids[k])
      local[[k]] <- igraph::layout_with_fr(subs[[k]])
    }
  }

  coords <- matrix(0, n, 2)
  if (identical(layout, "components")) {
    set.seed(seed)
    coords[unlist(fams), ] <- igraph::merge_coords(subs, local)
  } else {
    # Shelf packing: square cells of side sqrt(size) + 0.8, filled row by row.
    sizes <- lengths(fams)
    side <- sqrt(sizes) + 0.8
    width <- sqrt(sum(side^2)) * 1.1
    x <- 0
    y <- 0
    row_h <- 0
    for (k in order(-sizes)) {
      s <- side[k]
      if (x + s > width && x > 0) {
        x <- 0
        y <- y - row_h
        row_h <- 0
      }
      loc <- .rkn_centre_scale(local[[k]], (s - 0.8) / 2)
      coords[fams[[k]], ] <- cbind(loc[, 1] + x + s / 2, loc[, 2] + y - s / 2)
      x <- x + s
      row_h <- max(row_h, s)
    }
  }
  .rkn_centre_scale(coords, 1)
}

#' Internal helper: centre coordinates on the origin and scale them into a square
#'
#' Translates the midrange of each axis to zero and divides both axes by the same
#' factor, so the coordinates lie between -r and r and keep their aspect ratio. A single
#' point is returned at the origin.
#'
#' @keywords internal
#' @noRd
.rkn_centre_scale <- function(xy, r) {
  xy[, 1] <- xy[, 1] - (max(xy[, 1]) + min(xy[, 1])) / 2
  xy[, 2] <- xy[, 2] - (max(xy[, 2]) + min(xy[, 2])) / 2
  m <- max(abs(xy))
  if (m > 0) xy <- xy / m * r
  xy
}
