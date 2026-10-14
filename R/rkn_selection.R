#' Internal constructor for rkn_selection objects
#'
#' @keywords internal
#' @noRd
.new_rkn_selection <- function(kept, discarded, family_sizes, method, criteria, seed,
                               node_diagnostics, family_summary, graph, call) {
  structure(
    list(
      kept = kept,
      discarded = discarded,
      family_sizes = family_sizes,
      method = method,
      criteria = criteria,
      seed = seed,
      diagnostics = !is.null(node_diagnostics),
      node_diagnostics = node_diagnostics,
      family_summary = family_summary,
      graph = graph,
      call = call
    ),
    class = "rkn_selection"
  )
}

#' Print an rkn_selection object
#'
#' Prints a brief summary of a selection of unrelated individuals: how the selection
#' was made, how many individuals were retained, and how the network decomposes into
#' families.
#'
#' @param x An \code{rkn_selection} object, as returned by \code{\link{rkn_select_random}}
#'   or \code{\link{rkn_select_criteria}}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.rkn_selection <- function(x, ...) {
  n_ind <- length(x$kept) + length(x$discarded)
  sizes <- x$family_sizes

  header <- "Random selection"
  if (identical(x$method, "criteria")) {
    crit <- x$criteria$primary$type
    if (!is.null(x$criteria$secondary)) crit <- paste0(crit, ", then ", x$criteria$secondary$type)
    header <- sprintf("Criteria-based selection (%s)", crit)
  }
  if (!is.null(x$seed)) header <- paste0(header, ", seed ", x$seed)

  cat(sprintf("<rkn_selection> %s\n", header))
  cat(sprintf("Call: %s\n", paste(deparse(x$call), collapse = " ")))
  cat(sprintf("Individuals: %d | selected: %d (%.1f%%)\n",
              n_ind, length(x$kept), 100 * length(x$kept) / n_ind))
  cat(sprintf("Families: %d | size min/mean/max: %d / %.2f / %d\n",
              length(sizes), min(sizes), mean(sizes), max(sizes)))
  cat(sprintf("Families with up to 5 members: %d | with 6 or more: %d\n",
              sum(sizes <= 5), sum(sizes >= 6)))

  if (x$diagnostics) {
    n_mis <- sum(x$family_summary$used_largest_ivs %in% TRUE)
    cat(sprintf("Families resolved by enumerating maximum independent sets: %d\n", n_mis))
  } else {
    message("Node-level diagnostics were not computed (diagnostics = FALSE); ",
            "rerun with diagnostics = TRUE for a richer output.")
  }
  invisible(x)
}
