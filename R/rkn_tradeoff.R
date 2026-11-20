#' Compute the trade-off between selection size and composition
#'
#' @description
#' For a target composition of the whole selection, as in
#' \code{\link{rkn_select_global}}, shows what each individual given up allows: for each
#' number of individuals below the largest possible selection, the range of proportions
#' of the category, or of means of a numeric covariate, reached by the selections of that
#' size, and the one closest to the target. It helps to choose the \code{max_loss} of
#' \code{rkn_select_global()} before selecting.
#'
#' @inheritParams rkn_select_global
#' @param max_loss Either \code{NULL} (default) or a non-negative integer, the largest
#'   number of individuals given up covered by the table. With \code{NULL}, the table
#'   stops at the first number of individuals given up at which the target comes within
#'   reach (see Details), and at the latest at 5\% of the size of the largest selection,
#'   rounded up. With a number, it covers every number of individuals given up from 0 to
#'   \code{max_loss}, and at most the size of the largest selection.
#' @param seed Optional integer. Random numbers are used only in the families beyond the
#'   computational budget (see Details); a seed makes the table reproducible. If
#'   \code{NULL} (default), the current random number stream is used.
#'
#' @details
#' The largest possible selection keeps one maximum set of unrelated individuals per
#' family, and a family usually has several, with different numbers of members in the
#' category or different values of the covariate. The selections of the largest size
#' therefore reach a range of proportions or means, from the lowest to the highest. This
#' range is not an uncertainty: every proportion or mean in it is that of selections
#' that exist. Giving up individuals widens it, since members with any value can be left
#' out, and since smaller sets of unrelated members of a family can include individuals
#' that none of its maximum sets includes.
#'
#' For each number of individuals given up, the table reports the proportion or mean
#' closest to the target: the highest of the range when the target lies above it, the
#' lowest when the target lies below it, and, when the target lies within it, the
#' closest of those reached. A proportion among \eqn{n} individuals with a known value is
#' a multiple of \eqn{1/n}, so that it usually lies within \eqn{1/(2n)} of a target within
#' the range (the range may have gaps when few families offer alternatives). When two
#' values on either side of the target are equally close, the table reports the lower,
#' while \code{rkn_select_global()} draws among the selections reaching either. The
#' target comes within reach at the first number of individuals given up at which it
#' lies within the range.
#'
#' The distance from the target at the largest size is a property of the data, not an
#' approximation: it depends on the composition of the network, on the individuals
#' without relatives, who are always selected, and on the alternatives that each family
#' offers. \code{rkn_select_global(g, covariate, target, max_loss = L)} selects, among
#' the rows from 0 to \code{L}, the smallest distance and, among the rows reaching it,
#' the largest size. With a target beyond the range, such as \code{c(F = 1)}, it selects the
#' highest proportion of the category.
#'
#' For a numeric covariate, the means are those of the values rounded to multiples of
#' \code{resolution}, as in the search of \code{rkn_select_global()}: each lies within
#' half the resolution of the mean of the recorded values of the same selections, and
#' they are exact when no value is rounded, as for ages in whole years with the default
#' resolution.
#'
#' The table comes from the exact counts of \code{rkn_select_global()} (see its
#' Details), without drawing a selection. With \code{max_loss = NULL}, the table is
#' computed with no individual given up, then for 1, 2, 4 and so on, doubling until the
#' target comes within reach or the limit of 5\% is reached; the time grows faster than
#' the number of individuals given up, and much faster for a numeric covariate. The
#' table also stops at the last number computed if the next would exceed the memory
#' limit of \code{rkn_select_global()}, or its computational budget, of about half a
#' minute on a recent laptop; \code{print()} reports why. With a number as
#' \code{max_loss}, the function stops with an error beyond the memory limit, and a
#' message warns beyond the computational budget that the table may take several
#' minutes.
#'
#' The families beyond the computational budget of \code{rkn_select_global()} contribute
#' only maximum sets, and cannot give up individuals in the table: the range reachable
#' may then be wider than reported, and the number of selections is not available. Their
#' options come in part from random draws, so that \code{seed} makes the table
#' reproducible; \code{print()} reports their number.
#'
#' @return A data frame of class \code{rkn_tradeoff}, with one row per number of
#'   individuals given up and the columns:
#'   \describe{
#'     \item{loss}{Number of individuals given up, below the largest selection.}
#'     \item{size}{Size of the selection.}
#'     \item{proportion}{Proportion of the category closest to the target among the
#'       selections of that size (the lower of two equally close); \code{mean} for a
#'       numeric covariate.}
#'     \item{distance}{Absolute difference between \code{proportion} or \code{mean} and
#'       the target.}
#'     \item{min_proportion, max_proportion}{Lowest and highest proportions of the
#'       category among the selections of that size; \code{min_mean} and
#'       \code{max_mean} for a numeric covariate.}
#'     \item{log_count}{Natural logarithm of the number of selections of that size at
#'       the smallest distance from the target, on either side; \code{NA} when some
#'       family is beyond the computational budget.}
#'   }
#'   The proportions and means are computed among the selected individuals with a known
#'   value of the covariate, and are \code{NA} for a size at which no selection includes
#'   one. The attribute \code{"tradeoff"} is a list describing the table: \code{type}
#'   (\code{"proportion"} or \code{"mean"}); \code{covariate}; \code{category}, for a
#'   categorical covariate; \code{target} and \code{cohort}, whether the target is the
#'   composition of the network; \code{resolution}, \code{rounded} and
#'   \code{coarsened}, for a numeric covariate, whether some value was rounded and
#'   whether the default resolution was coarsened to fit the memory budget;
#'   \code{max_size}, the size of the largest selection; \code{max_loss}, as given; \code{reached}, the number of individuals
#'   given up at which the target comes within reach (\code{NA} if it does not within the
#'   table); \code{limited}, why the default table stopped before the target came within
#'   reach, \code{"memory"} or \code{"work"} for the memory limit or the computational
#'   budget (\code{NA} otherwise); \code{n_missing}, the number of individuals of the
#'   network with a missing value; and \code{n_inexact}, the number of families beyond
#'   the computational budget.
#'
#'   \code{print()} shows the first rows and those around \code{reached}, with the
#'   number of selections as a count. \code{plot()} draws, against the number of
#'   individuals given up, the range of proportions or means (shaded), the one closest to
#'   the target (line and points) and the target (dashed), with the size of the
#'   selection on the top axis and, if the target comes within reach after giving up
#'   individuals, a dotted line at that number.
#'
#' @seealso \code{\link{rkn_select_global}}.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' # How many individuals to give up for 70% women
#' tr <- rkn_tradeoff(g, "sex", c(F = 0.7))
#' tr
#' plot(tr)
#'
#' # The selection with that composition
#' rkn_select_global(g, "sex", c(F = 0.7), max_loss = 1, seed = 1)
#'
#' # Every number of individuals given up, from 0 to 4
#' rkn_tradeoff(g, "sex", c(F = 0.7), max_loss = 4)
#'
#' # The proportion of women of the network is within reach at the largest size
#' rkn_tradeoff(g, "sex")
#'
#' # A mean age of 45 years
#' plot(rkn_tradeoff(g, "age", 45, max_loss = 4))
#'
#' @export
rkn_tradeoff <- function(g, covariate, target = NULL, max_loss = NULL, resolution = NULL, seed = NULL) {

  if (missing(covariate)) covariate <- NULL
  .rkn_global_check_input(g, covariate)
  if (!is.null(max_loss) && (!is.numeric(max_loss) || length(max_loss) != 1 || is.na(max_loss) ||
                             max_loss < 0 || max_loss != round(max_loss))) {
    stop("`max_loss` must be NULL, to extend the table until the target comes within reach (up to 5% ",
         "of the largest selection), or a single non-negative integer.", call. = FALSE)
  }
  .rkn_check_seed(seed)
  tgt <- .rkn_global_target(igraph::vertex_attr(g, covariate), target, resolution)
  if (!is.null(seed)) {
    withr::local_preserve_seed()
    set.seed(seed)
  }
  p <- tgt$p
  within_reach <- function(rows) {
    which(rows$min_value <= p + 1e-12 & rows$max_value >= p - 1e-12)
  }

  # With max_loss = NULL the table starts with no individual given up and is extended
  # by doubling their number; the families whose options do not change are reused.
  k <- if (is.null(max_loss)) 0 else max_loss
  fam <- NULL
  rows <- NULL
  used <- NULL
  limited <- NA_character_
  coarsened <- FALSE
  repeat {
    fam <- .rkn_global_families(g, tgt, k, fam)
    limit <- min(if (is.null(max_loss)) ceiling(0.05 * fam$max_size) else max_loss, fam$max_size)
    # Beyond the computational budget, the first table is computed with a message,
    # and the extension of the default table stops. Beyond the memory budget, the default
    # resolution of a numeric covariate is coarsened for the first table, as in
    # rkn_select_global(), and the extension of the default table stops.
    res <- tryCatch(.rkn_global_select(fam$fams, min(k, limit), draw = FALSE,
                                       long = if (is.null(max_loss) && !is.null(rows)) "stop" else "message"),
                    rkn_table_limit = function(e) {
                      if (!is.null(rows)) return("memory")
                      if (is.null(resolution) && tgt$type == "mean") e else stop(e)
                    },
                    rkn_work_limit = function(e) "work")
    if (inherits(res, "rkn_table_limit")) {
      tgt <- .rkn_global_target(igraph::vertex_attr(g, covariate), target,
                                .rkn_coarser_resolution(tgt$resolution, res$need / res$max_cells))
      fam <- NULL
      coarsened <- TRUE
      next
    }
    if (is.character(res)) {
      limited <- res
      break
    }
    rows <- .rkn_tradeoff_rows(res, tgt)
    used <- fam
    if (!is.null(max_loss) || length(within_reach(rows)) > 0 || k >= limit) break
    k <- min(max(1, 2 * k), limit)
  }

  if (all(is.na(rows$value))) {
    last <- max(rows$loss)
    stop(sprintf("No selection giving up at most %d %s includes an individual with a known value of \"%s\".",
                 last, if (last == 1) "individual" else "individuals", covariate), call. = FALSE)
  }
  within <- within_reach(rows)
  reached <- if (length(within) > 0) rows$loss[within[1]] else NA_integer_
  if (is.null(max_loss) && !is.na(reached)) rows <- rows[seq_len(reached + 1), ]
  names(rows)[3:6] <- .rkn_tradeoff_columns(tgt$type)
  attr(rows, "tradeoff") <- list(
    type = tgt$type, covariate = covariate, category = tgt$category, target = p, cohort = tgt$cohort,
    resolution = if (tgt$type == "mean") tgt$resolution, rounded = isTRUE(tgt$rounded), coarsened = coarsened,
    max_size = used$max_size, max_loss = max_loss, reached = reached, limited = limited,
    n_missing = sum(tgt$miss), n_inexact = sum(!used$exact)
  )
  class(rows) <- c("rkn_tradeoff", "data.frame")
  rows
}

#' Internal helper: the rows of a trade-off table
#'
#' @param res The result of \code{.rkn_global_select(draw = FALSE)}.
#' @param tgt The result of \code{.rkn_global_target()}.
#' @return A data frame with one row per number of individuals given up, with the
#'   columns described in \code{rkn_tradeoff()}, the proportion or mean and its range
#'   being named \code{value}, \code{min_value} and \code{max_value}.
#'
#' @keywords internal
#' @noRd
.rkn_tradeoff_rows <- function(res, tgt) {
  tab <- res$table
  loss <- seq_len(dim(tab)[1]) - 1L
  distance <- .rkn_global_distance(tgt)
  out <- data.frame(loss = loss, size = as.integer(res$alpha - loss), value = NA_real_,
                    distance = NA_real_, min_value = NA_real_, max_value = NA_real_, log_count = NA_real_)
  for (j in loss) {
    cells <- which(tab[j + 1, , , drop = FALSE] > -Inf, arr.ind = TRUE)
    M <- res$m_base + cells[, 3] - 1
    known <- res$alpha - j - M
    cells <- cells[known > 0, , drop = FALSE]
    M <- M[known > 0]
    known <- known[known > 0]
    if (length(known) == 0) next
    C <- res$c_base + cells[, 2] - 1
    value <- tgt$lowest + tgt$resolution * C / known
    d <- distance(C, M, res$alpha - j)
    best <- d <= min(d) + 1e-12
    out$value[j + 1] <- min(value[best])
    out$distance[j + 1] <- min(d)
    out$min_value[j + 1] <- min(value)
    out$max_value[j + 1] <- max(value)
    if (res$exact) {
      lw <- tab[cbind(j + 1, cells[best, 2], cells[best, 3])]
      out$log_count[j + 1] <- max(lw) + log(sum(exp(lw - max(lw))))
    }
  }
  out
}

#' Internal helper: names of the columns of a trade-off table
#'
#' @param type \code{"proportion"} or \code{"mean"}.
#' @return The names of the closest value, the distance and the range.
#'
#' @keywords internal
#' @noRd
.rkn_tradeoff_columns <- function(type) {
  c(type, "distance", paste0(c("min_", "max_"), type))
}

#' @param x An \code{rkn_tradeoff} object.
#' @param ... Not used.
#'
#' @rdname rkn_tradeoff
#' @export
print.rkn_tradeoff <- function(x, ...) {
  info <- attr(x, "tradeoff")
  if (is.null(info)) return(NextMethod())
  cols <- .rkn_tradeoff_columns(info$type)
  proportion <- info$type == "proportion"
  fmt <- function(v) {
    out <- if (proportion) sprintf("%.4f", v) else .rkn_format_mean(v, info$resolution)
    out[is.na(v)] <- "NA"
    out
  }
  label <- if (proportion) sprintf("Proportion of %s = %s", info$covariate, info$category) else
    sprintf("Mean of %s", info$covariate)
  cat(sprintf("<rkn_tradeoff> %s by individuals given up (target %s%s)\n", label,
              if (proportion) sprintf("%.4g", info$target) else format(info$target, digits = 6),
              if (!info$cohort) "" else if (proportion) ", the proportion of the network" else ", the mean of the network"))
  first <- which(x$loss == 0)
  cat(sprintf("Largest selection: %d individuals%s\n", info$max_size,
              if (length(first) == 1) sprintf(", %s from %s to %s",
                                               if (proportion) paste("proportion of", info$category) else
                                                 paste("mean of", info$covariate),
                                               fmt(x[[cols[3]]][first]), fmt(x[[cols[4]]][first])) else ""))
  last <- x$loss[nrow(x)]
  share <- function(n) sprintf("%.1f%% of the largest selection", 100 * n / info$max_size)
  if (!is.na(info$reached) && info$reached == 0) {
    cat("Target within reach at the largest selection: no individual needs to be given up\n")
  } else if (!is.na(info$reached)) {
    cat(sprintf("Target within reach by giving up %d (%s), a selection of %d\n",
                info$reached, share(info$reached), info$max_size - info$reached))
  } else if (!is.na(info$limited)) {
    cat(sprintf("Target not within reach by giving up up to %d (%s); a longer table would exceed %s\n",
                last, share(last),
                if (info$limited == "memory") "the memory limit" else
                  "the computational budget: set max_loss to extend it"))
  } else if (is.null(info$max_loss)) {
    cat(sprintf(paste0("Target not within reach by giving up up to %d (5%% of the largest selection, rounded up); ",
                       "set max_loss to extend the table\n"), last))
  } else {
    cat(sprintf("Target not within reach by giving up up to %d (%s)\n", last, share(last)))
  }
  if (info$rounded) {
    cat(sprintf("Values of %s rounded to multiples of %s%s: each mean is within %s of that of the recorded values\n",
                info$covariate, format(info$resolution),
                if (isTRUE(info$coarsened)) ", coarser than the default to fit the memory budget" else "",
                format(info$resolution / 2)))
  }
  if (info$n_missing > 0) {
    cat(sprintf("%s among the selected individuals with a known value of %s\n",
                if (proportion) "Proportions" else "Means", info$covariate))
  }
  if (info$n_inexact > 0) {
    cat(sprintf(paste0("Families beyond the computational budget: %d; they cannot give up individuals here, ",
                       "the %s reachable may be wider, and the number of selections is not available\n"),
                info$n_inexact, if (proportion) "proportions" else "means"))
  }

  n <- nrow(x)
  shown <- seq_len(n)
  if (n > 20) {
    around <- if (is.na(info$reached)) integer(0) else info$reached + 0:2
    shown <- sort(unique(c(1:6, around[around <= n], n)))
  }
  selections <- rep("NA", n)
  for (i in which(!is.na(x$log_count))) selections[i] <- .rkn_format_log10(x$log_count[i] / log(10))
  out <- data.frame(loss = x$loss, size = x$size, fmt(x[[cols[1]]]), fmt(x$distance), fmt(x[[cols[3]]]),
                    fmt(x[[cols[4]]]), selections)
  names(out)[3:7] <- c(cols, "selections")
  print(out[shown, ], row.names = FALSE)
  if (length(shown) < n) {
    cat(sprintf("(%d of %d rows; the object is a data frame with all of them)\n", length(shown), n))
  }
  invisible(x)
}

#' @param main Title of the plot; \code{NULL} (default) describes the covariate, the
#'   category if any, and the target.
#'
#' @rdname rkn_tradeoff
#' @export
plot.rkn_tradeoff <- function(x, main = NULL, ...) {
  info <- attr(x, "tradeoff")
  if (is.null(info)) stop("`x` has lost the attributes of the object returned by rkn_tradeoff().", call. = FALSE)
  cols <- .rkn_tradeoff_columns(info$type)
  ok <- !is.na(x[[cols[1]]])
  loss <- x$loss[ok]
  closest <- x[[cols[1]]][ok]
  lo <- x[[cols[3]]][ok]
  hi <- x[[cols[4]]][ok]
  p <- info$target
  label <- if (info$type == "proportion") sprintf("Proportion of %s = %s", info$covariate, info$category) else
    sprintf("Mean of %s", info$covariate)
  if (is.null(main)) main <- sprintf("%s by individuals given up (target %s)", label, format(p, digits = 4))
  single <- length(loss) == 1
  xlim <- if (single) loss + c(-1, 1) else range(loss)
  ylim <- range(c(lo, hi, p))

  op <- graphics::par(mar = c(6.2, 4.2, 4.2, 1), mgp = c(2.4, 0.7, 0), las = 1)
  on.exit(graphics::par(op))
  graphics::plot(loss, closest, type = "n", xlim = xlim, ylim = ylim, xaxt = "n", bty = "l",
                 xlab = "Individuals given up", ylab = label)
  at <- if (single) loss else unique(round(pretty(loss)))
  at <- at[at >= min(loss) & at <= max(loss)]
  graphics::axis(1, at = at)
  graphics::axis(3, at = at, labels = info$max_size - at)
  graphics::mtext("Selection size", side = 3, line = 1.6, cex = 0.8)
  graphics::title(main, line = 2.6, cex.main = 1)
  if (single) {
    graphics::rect(loss - 0.15, lo, loss + 0.15, hi, col = "#C6DBEF", border = NA)
  } else {
    graphics::polygon(c(loss, rev(loss)), c(lo, rev(hi)), col = "#C6DBEF", border = NA)
  }
  graphics::abline(h = p, lty = 2, col = "#E69F00", lwd = 1.6)
  graphics::lines(loss, closest, col = "#0072B2", lwd = 1.8)
  graphics::points(loss, closest, pch = 19, col = "#0072B2", cex = if (length(loss) > 30) 0.45 else 1)
  if (!is.na(info$reached) && info$reached > 0) {
    r <- info$reached
    right <- r > mean(xlim)
    graphics::abline(v = r, lty = 3, col = "grey40")
    graphics::text(r, ylim[1] + 0.06 * diff(ylim),
                   sprintf(if (right) "target within reach by giving up %d " else " target within reach by giving up %d", r),
                   adj = if (right) 1 else 0, cex = 0.75, col = "grey25")
  }

  # The legend sits below the axis title, at a fixed number of lines of the margin.
  usr <- graphics::par("usr")
  inch <- graphics::par("mai")[1] / graphics::par("mar")[1]
  y <- usr[3] - 4.4 * inch * diff(usr[3:4]) / graphics::par("pin")[2]
  graphics::legend(mean(usr[1:2]), y, xjust = 0.5, yjust = 0.5, xpd = NA, horiz = TRUE, bty = "n", cex = 0.8,
                   legend = c(if (info$type == "proportion") "proportions reachable" else "means reachable",
                              "closest to the target", "target"),
                   fill = c("#C6DBEF", NA, NA), border = NA, lty = c(NA, 1, 2), lwd = c(NA, 1.8, 1.6),
                   pch = c(NA, 19, NA), col = c(NA, "#0072B2", "#E69F00"))
  invisible(x)
}
