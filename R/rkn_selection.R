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

#' Internal helper: one-line description of how a selection was made
#'
#' @keywords internal
#' @noRd
.rkn_selection_header <- function(x) {
  header <- "Random selection"
  if (identical(x$method, "criteria")) {
    crit <- x$criteria$primary$type
    if (!is.null(x$criteria$secondary)) crit <- paste0(crit, ", then ", x$criteria$secondary$type)
    header <- sprintf("Criteria-based selection (%s)", crit)
  }
  if (!is.null(x$seed)) header <- paste0(header, ", seed ", x$seed)
  header
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

  cat(sprintf("<rkn_selection> %s\n", .rkn_selection_header(x)))
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

#' Summarise an rkn_selection object
#'
#' Computes summary tables of a selection of unrelated individuals: selection counts by
#' family size and, when diagnostics were computed, the structures recognised in the
#' families and how each individual was decided.
#'
#' @param object An \code{rkn_selection} object, as returned by
#'   \code{\link{rkn_select_random}} or \code{\link{rkn_select_criteria}}.
#' @param ... Not used.
#'
#' @return An object of class \code{summary.rkn_selection}: a list with elements
#' \describe{
#'   \item{\code{header}}{a one-line description of how the selection was made;}
#'   \item{\code{overall}}{a one-row data frame with the numbers of individuals, of
#'     selected individuals and of families;}
#'   \item{\code{by_size}}{a data frame with one row per family size, giving the number
#'     of families, of individuals and of selected individuals;}
#'   \item{\code{structures}}{a data frame with the number of families per recognised
#'     structure, or \code{NULL} if diagnostics were not computed;}
#'   \item{\code{decisions}}{a contingency table of individuals by type of decision
#'     (a priori or among candidates) and outcome (kept or discarded), or \code{NULL}
#'     if diagnostics were not computed.}
#' }
#'
#' @export
summary.rkn_selection <- function(object, ...) {
  sizes <- object$family_sizes
  membership <- igraph::components(object$graph)$membership
  names(membership) <- igraph::V(object$graph)$name
  kept_per_family <- tabulate(membership[object$kept], nbins = length(sizes))

  size_levels <- sort(unique(sizes))
  n_families <- integer(length(size_levels))
  n_selected <- integer(length(size_levels))
  for (k in seq_along(size_levels)) {
    in_level <- sizes == size_levels[k]
    n_families[k] <- sum(in_level)
    n_selected[k] <- sum(kept_per_family[in_level])
  }
  by_size <- .rkn_df(
    family_size = size_levels,
    n_families = n_families,
    n_individuals = size_levels * n_families,
    n_selected = n_selected
  )

  n_ind <- length(object$kept) + length(object$discarded)
  overall <- .rkn_df(
    n_individuals = n_ind,
    n_selected = length(object$kept),
    n_families = length(sizes)
  )

  structures <- NULL
  decisions <- NULL
  if (object$diagnostics) {
    counts <- sort(table(object$family_summary$structure_name), decreasing = TRUE)
    structures <- .rkn_df(structure = names(counts), n_families = as.integer(counts))
    decisions <- table(decision = object$node_diagnostics$pre_class,
                       outcome = object$node_diagnostics$post_class)
  }

  structure(
    list(header = .rkn_selection_header(object), overall = overall, by_size = by_size,
         structures = structures, decisions = decisions),
    class = "summary.rkn_selection"
  )
}

#' Print a summary of an rkn_selection object
#'
#' @param x A \code{summary.rkn_selection} object, as returned by
#'   \code{\link{summary.rkn_selection}}.
#' @param max_structures Maximum number of recognised structures to list.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @method print summary.rkn_selection
#' @export
print.summary.rkn_selection <- function(x, max_structures = 10, ...) {
  o <- x$overall
  cat(sprintf("<rkn_selection summary> %s\n", x$header))
  cat(sprintf("Individuals: %d | selected: %d (%.1f%%) | families: %d\n\n",
              o$n_individuals, o$n_selected, 100 * o$n_selected / o$n_individuals, o$n_families))

  cat("Selection by family size:\n")
  print(x$by_size, row.names = FALSE)

  if (is.null(x$structures)) {
    message("Node-level diagnostics were not computed (diagnostics = FALSE); ",
            "rerun with diagnostics = TRUE for a richer output.")
    return(invisible(x))
  }

  n_struct <- nrow(x$structures)
  cat(sprintf("\nFamilies by recognised structure (%s):\n",
              if (n_struct > max_structures) sprintf("top %d of %d", max_structures, n_struct) else "all"))
  print(x$structures[seq_len(min(n_struct, max_structures)), , drop = FALSE], row.names = FALSE)

  cat("\nIndividuals by type of decision and outcome:\n")
  print(x$decisions)
  invisible(x)
}

#' Plot an rkn_selection object
#'
#' Draws the relatedness network of a selection, with the vertices coloured according
#' to \code{type}.
#'
#' @param x An \code{rkn_selection} object, as returned by \code{\link{rkn_select_random}}
#'   or \code{\link{rkn_select_criteria}}.
#' @param type One of \code{"post"} (default; selected vs discarded individuals, see
#'   \code{\link{rkn_plot_post}}), \code{"diagnostics"} (individuals coloured by type of
#'   decision, see \code{\link{rkn_plot_diagnostics}}; requires diagnostics) or
#'   \code{"pre"} (families before the selection, see \code{\link{rkn_plot_pre}}).
#' @param ... Further arguments passed to the plotting function selected by
#'   \code{type}, for instance \code{layout} or \code{vertex_size}.
#'
#' @return \code{x}, invisibly.
#'
#' @export
plot.rkn_selection <- function(x, type = c("post", "diagnostics", "pre"), ...) {
  type <- match.arg(type)
  if (identical(type, "diagnostics") && !x$diagnostics) {
    message("Node-level diagnostics were not computed (diagnostics = FALSE); plotting ",
            "selected vs discarded instead. Rerun with diagnostics = TRUE for the diagnostics plot.")
    type <- "post"
  }

  if (identical(type, "post")) {
    rkn_plot_post(x$graph, kept = x$kept, ...)
  } else if (identical(type, "diagnostics")) {
    rkn_plot_diagnostics(x$graph, res = x, ...)
  } else {
    rkn_plot_pre(x$graph, ...)
  }
  invisible(x)
}
