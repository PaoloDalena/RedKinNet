#' Selection with a target composition of the whole selection
#'
#' @description
#' Selects a maximum set of unrelated individuals whose composition, over the whole
#' selection, is as close as possible to a target: the proportion of one category of a
#' categorical covariate, such as the proportion of women, or the mean of a numeric
#' covariate, such as the mean age. Unlike \code{\link{rkn_select_criteria}}, which
#' applies its criteria within each family, the families are chosen jointly: the choice
#' in one family depends on what the others contribute. For instance, with ten unrelated
#' women and ten related pairs of a woman and a man, choosing within each pair gives
#' fifteen women on average, while choosing the man in every pair gives ten women and ten
#' men.
#'
#' @inheritParams rkn_select_random
#' @param covariate Character. Name of a vertex attribute of the network, either
#'   categorical (character, factor or logical), such as \code{"sex"} for
#'   \code{V(g)$sex}, or numeric, such as \code{"age"}.
#' @param target For a categorical covariate, the target proportion of one category
#'   among the selected individuals, named by the category, as \code{c(F = 0.5)}; for a
#'   covariate with two categories, the proportions of both may be given, summing to 1,
#'   and the first is used. For a numeric covariate, the target mean, a single unnamed
#'   number, as \code{45}. If \code{NULL} (default), the target is the proportion of the
#'   first category, in alphabetical order, or the mean, in the whole network, so that
#'   the selection keeps the composition of the cohort; for a categorical covariate this
#'   requires two categories.
#' @param max_loss Non-negative integer. Largest number of individuals that may be
#'   selected below the maximum, to bring the composition closer to the target. With
#'   \code{0} (default) the selection keeps the maximum size. \code{\link{rkn_tradeoff}}
#'   shows, before selecting, the proportion or mean that each number of individuals
#'   given up allows.
#' @param resolution For numeric covariates only, a positive number: the values are
#'   rounded to multiples of \code{resolution} above the smallest value to search the
#'   selections (see Details). If \code{NULL} (default), integer values spanning at most
#'   100 are not rounded (a resolution of 1, as for ages in years or scores); other values
#'   are rounded to the smallest of 1, 2 or 5 times a power of 10 that divides their range
#'   into at most 100 steps, as 0.5 for values from 16 to 45, or 2000 for values from 0
#'   to 200,000.
#'
#' @details
#' The selection is chosen in three steps: first the proportion or mean closest to the
#' target; then, among the selections that reach it, the largest; then one of the
#' selections with that composition and that size, drawn uniformly at random. With
#' \code{max_loss > 0}, individuals are given up even for a small improvement of the
#' composition. The proportion or the mean is computed among the selected individuals
#' with a known value of the covariate; the number of selected individuals with a
#' missing value is reported by \code{print()}.
#'
#' The optimisation is exact. For each family, the exact algorithm described in the
#' section Search over maximum sets of \code{\link{rkn_selection}} counts its sets of
#' unrelated members by size and by number of members in the category, or by sum of the
#' values; dynamic programming then combines the families. The counts are kept on the
#' logarithmic scale, so that rare combinations, as for an extreme target, are not lost,
#' and \code{print()} reports the number of selections reaching the best composition.
#'
#' For a numeric covariate, counting the selections by sum of the values requires values
#' on a grid: with values recorded to the day, almost every selection would have its own
#' sum. The values are therefore rounded to multiples of \code{resolution} for the search
#' only, and \code{print()} reports the mean of the recorded values of the selected
#' individuals. Since the rounding moves the mean of any selection by at most half the
#' resolution, the distance of the selected mean from the target exceeds the smallest
#' possible by at most \code{resolution}; in practice the rounding errors largely cancel,
#' and the mean of \eqn{n} individuals moves by about
#' \eqn{\mathrm{resolution}/\sqrt{12 n}}{resolution / sqrt(12 n)}. The time grows with
#' the number of steps over the range of the values, and quickly with \code{max_loss}: a
#' coarser resolution shortens it.
#'
#' The computation of a family is not attempted when it would exceed a budget of memory
#' and time, as for the largest families of genetic isolates at the third degree of
#' relationship. Such a family contributes the maximum sets with the smallest and the
#' largest number of members in the category, or sum of the values, found exactly, and
#' those of 100 uniform draws; a family too intricate even for this is selected as by
#' \code{\link{rkn_select_random}}. The selection then keeps the maximum size, but the
#' composition may not be the best possible and the draw is not uniform within these
#' families, whose number \code{print()} reports. If the tables that combine the
#' families would exceed about 80 MB, which requires many individuals with a missing
#' value, a large \code{max_loss} or, for a numeric covariate, a fine resolution, the
#' selection stops with an error before allocating them. The work of these tables is
#' also estimated beforehand: beyond about half a minute on a recent laptop, a message
#' warns that the selection may take several minutes, and a smaller \code{max_loss} or a
#' coarser resolution shortens it.
#'
#' @return An object of class \code{\link{rkn_selection}} with \code{method = "global"};
#'   its element \code{global} describes the composition reached (see
#'   \code{\link{rkn_selection}}).
#'
#' @seealso \code{\link{rkn_tradeoff}} for the proportions or means reachable by each
#'   number of individuals given up; \code{\link{rkn_selection}},
#'   \code{\link{rkn_select_criteria}}, \code{\link{rkn_select_random}}.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' # Half women among the selected
#' sel <- rkn_select_global(g, "sex", c(F = 0.5), seed = 1)
#' sel
#' sel$kept
#'
#' # The proportion of women of the whole network
#' rkn_select_global(g, "sex", seed = 1)
#'
#' # Give up to two individuals to get closer to 70% women
#' rkn_select_global(g, "sex", c(F = 0.7), max_loss = 2, seed = 1)
#'
#' # A mean age of 45 years
#' rkn_select_global(g, "age", 45, seed = 1)
#'
#' @export
rkn_select_global <- function(g, covariate, target = NULL, max_loss = 0, resolution = NULL, seed = NULL,
                              diagnostics = TRUE) {

  if (missing(covariate)) covariate <- NULL
  .rkn_global_check_input(g, covariate)
  if (!is.numeric(max_loss) || length(max_loss) != 1 || is.na(max_loss) || max_loss < 0 ||
      max_loss != round(max_loss)) {
    stop("`max_loss` must be a single non-negative integer.")
  }
  values <- igraph::vertex_attr(g, covariate)
  tgt <- .rkn_global_target(values, target, resolution)

  if (!is.null(seed)) {
    withr::local_preserve_seed()
    set.seed(seed)
  }

  ids <- igraph::V(g)$name
  fam <- .rkn_global_families(g, tgt, max_loss)
  comp <- fam$comp

  # The network cannot give up more individuals than its maximum size.
  res <- .rkn_global_select(fam$fams, min(max_loss, fam$max_size), .rkn_global_distance(tgt))
  sets <- .rkn_global_realise(fam$fams, res$choice)
  selected <- rep(FALSE, length(ids))
  for (i in seq_along(sets)) selected[fam$members[[i]][sets[[i]]]] <- TRUE
  known <- selected & tgt$miss == 0

  global <- list(
    type = tgt$type, covariate = covariate, target = tgt$p, cohort = tgt$cohort,
    n_known = res$size - res$m, n_missing = res$m, size = res$size, max_size = res$size + res$j,
    max_loss = max_loss, log_count = res$log_count, exact = res$exact, n_inexact = sum(!fam$exact)
  )
  if (tgt$type == "proportion") {
    global$category <- tgt$category
    global$proportion <- res$c / (res$size - res$m)
    global$n_category <- res$c
  } else {
    global$mean <- mean(values[known])
    global$resolution <- tgt$resolution
    global$rounded <- tgt$rounded
  }

  node_diagnostics <- NULL
  family_summary <- NULL
  if (diagnostics) {
    family_id <- comp$membership
    family_size <- comp$csize[family_id]
    node_diagnostics <- .rkn_df(
      id = ids, pre_class = rep(NA_character_, length(ids)),
      post_class = ifelse(selected, "kept", "discarded"), is_selected = selected,
      pre_reason = rep(NA_character_, length(ids)), family_id = family_id, family_size = family_size
    )
    family_summary <- .rkn_df(
      family_id = seq_len(comp$no), family_size = comp$csize,
      selector = rep("global", comp$no), exact = fam$exact,
      n_kept = tabulate(family_id[selected], nbins = comp$no),
      n_discarded = tabulate(family_id[!selected], nbins = comp$no)
    )
    if (tgt$type == "proportion") {
      family_summary$n_category <- tabulate(family_id[selected & tgt$cat == 1], nbins = comp$no)
    } else {
      sums <- numeric(comp$no)
      for (i in which(known)) sums[family_id[i]] <- sums[family_id[i]] + values[i]
      n <- tabulate(family_id[known], nbins = comp$no)
      family_summary$mean <- ifelse(n > 0, sums / n, NA_real_)
    }
    family_summary$n_missing <- tabulate(family_id[selected & tgt$miss == 1], nbins = comp$no)
  }

  .new_rkn_selection(
    kept = ids[selected], discarded = ids[!selected], family_sizes = comp$csize,
    method = "global", criteria = NULL, seed = seed,
    node_diagnostics = node_diagnostics, family_summary = family_summary, graph = g,
    call = match.call(), global = global
  )
}

#' Internal helper: checks of the network and covariate of a global selection
#'
#' @param g The network.
#' @param covariate The name of the covariate, or \code{NULL} if not given.
#' @return \code{NULL}, invisibly; stops with an error on an invalid input.
#'
#' @keywords internal
#' @noRd
.rkn_global_check_input <- function(g, covariate) {
  if (!igraph::is_igraph(g)) stop("g must be an igraph object.")
  if (igraph::vcount(g) == 0) stop("g must contain at least one vertex.")
  if (is.null(igraph::V(g)$name)) stop("V(g)$name must be defined (vertex identifiers are required).")
  reserved <- intersect(igraph::vertex_attr_names(g), .rkn_reserved_attributes)
  if (length(reserved) > 0) stop(.rkn_reserved_message(reserved))
  if (!is.character(covariate) || length(covariate) != 1) {
    stop("`covariate` must be the name of a vertex attribute, as \"sex\".")
  }
  if (!(covariate %in% igraph::vertex_attr_names(g))) {
    stop(sprintf("The network has no vertex attribute \"%s\".", covariate))
  }
  invisible(NULL)
}

#' Internal helper: the options of the families of a global selection
#'
#' Groups the individuals without relatives by value and computes the options of every
#' family (\code{.rkn_global_options()}). A family cannot give up more members than it
#' has, so that its options are computed for at most that many: a larger \code{k} would
#' only enlarge its table. Nor can it give up more than the size of its maximum sets, so
#' that options computed for at least that many are final.
#'
#' @param g The network.
#' @param tgt The result of \code{.rkn_global_target()}.
#' @param k Largest number of individuals below the maximum size of the network.
#' @param previous A result of this function for another \code{k}: the options of the
#'   families for which the limit is unchanged, or already final, are reused.
#' @return A list with \code{fams} (the options, for the groups of individuals without
#'   relatives first and then for the families), \code{members} (their vertices),
#'   \code{comp} (the components of the network), \code{options} and \code{cap} (the
#'   options of each component and the limit they were computed for), \code{exact}
#'   (per component) and \code{max_size}, the size of the largest selection.
#'
#' @keywords internal
#' @noRd
.rkn_global_families <- function(g, tgt, k, previous = NULL) {
  comp <- if (is.null(previous)) igraph::components(g) else previous$comp
  single <- which(comp$csize[comp$membership] == 1)
  fams <- list()
  members <- list()
  for (grp in .rkn_global_singletons(tgt$cat[single], tgt$miss[single], k)) {
    fams[[length(fams) + 1]] <- grp
    members[[length(members) + 1]] <- single[grp$members]
  }
  cap <- pmin(k, comp$csize)
  options <- vector("list", comp$no)
  exact <- rep(TRUE, comp$no)
  for (f in which(comp$csize > 1)) {
    v <- which(comp$membership == f)
    if (!is.null(previous) && (previous$cap[f] == cap[f] || previous$cap[f] >= previous$options[[f]]$alpha)) {
      options[[f]] <- previous$options[[f]]
      cap[f] <- previous$cap[f]
    } else {
      g_sub <- igraph::induced_subgraph(g, v)
      fallback <- function() igraph::V(g_sub)$name %in% rkn_select_random(g_sub, diagnostics = FALSE)$kept
      options[[f]] <- .rkn_global_options(.rkn_adjacency(g_sub), tgt$cat[v], tgt$miss[v], cap[f],
                                          fallback = fallback)
    }
    fams[[length(fams) + 1]] <- options[[f]]
    members[[length(members) + 1]] <- v
    exact[f] <- options[[f]]$exact
  }
  max_size <- 0
  for (fa in fams) max_size <- max_size + fa$alpha
  list(fams = fams, members = members, comp = comp, options = options, cap = cap, exact = exact,
       max_size = max_size)
}

#' Internal helper: target of a global selection
#'
#' For a categorical covariate, the target is the proportion of a category; for a
#' numeric covariate, the target is the mean, and the values are rounded to multiples of
#' the resolution above the smallest value, as integer steps.
#'
#' @param x The covariate, one value per individual.
#' @param target The \code{target} argument of \code{rkn_select_global()}.
#' @param resolution The \code{resolution} argument of \code{rkn_select_global()}.
#' @return A list with \code{type} (\code{"proportion"} or \code{"mean"}), \code{p}
#'   (the target), \code{cohort} (whether the target is that of the network),
#'   \code{cat} (1 in the category, or the steps of the values), \code{miss} (1 for a
#'   missing value), and \code{lowest} and \code{resolution}, by which a mean of
#'   \code{cat} converts to the scale of the covariate (0 and 1 for a proportion);
#'   \code{category} for a categorical covariate; \code{rounded} (whether some value
#'   was rounded) for a numeric one.
#'
#' @keywords internal
#' @noRd
.rkn_global_target <- function(x, target, resolution = NULL) {
  if (is.numeric(x)) return(.rkn_global_target_mean(x, target, resolution))
  if (!is.null(resolution)) stop("`resolution` applies to numeric covariates only.")
  known <- !is.na(x)
  x <- as.character(x)
  categories <- sort(unique(x[known]))
  if (length(categories) == 0) stop("The covariate has no known value.")
  if (is.null(target)) {
    if (length(categories) > 2) {
      stop("The covariate has more than two categories: name one in `target` with its target ",
           "proportion, as c(A = 0.3).")
    }
    category <- categories[1]
    p <- mean(x[known] == category)
  } else {
    if (!is.numeric(target) || is.null(names(target)) || any(names(target) == "") || anyNA(target)) {
      stop("`target` must name the category with its target proportion, as c(F = 0.5).")
    }
    if (any(target < 0 | target > 1)) stop("The proportions in `target` must lie between 0 and 1.")
    if (length(target) > 2 || (length(target) == 2 && length(categories) > 2)) {
      stop("`target` sets the proportion of one category against the others: give one named ",
           "proportion, or two for a covariate with two categories.")
    }
    if (length(target) == 2 && abs(sum(target) - 1) > 1e-8) {
      stop("The two proportions in `target` must sum to 1.")
    }
    absent <- setdiff(names(target), categories)
    if (length(absent) > 0) {
      stop(sprintf("Category %s not found among the values of the covariate.",
                   paste0("\"", absent, "\"", collapse = ", ")))
    }
    category <- names(target)[1]
    p <- unname(target[1])
  }
  list(type = "proportion", category = category, p = p, cohort = is.null(target),
       cat = as.numeric(known & x == category), miss = as.numeric(!known), lowest = 0, resolution = 1)
}

#' Internal helper: target mean of a numeric covariate
#'
#' @inheritParams .rkn_global_target
#' @return As \code{.rkn_global_target()}.
#'
#' @keywords internal
#' @noRd
.rkn_global_target_mean <- function(x, target, resolution) {
  known <- !is.na(x)
  if (!any(known)) stop("The covariate has no known value.")
  if (any(is.infinite(x[known]))) stop("The covariate has infinite values.")
  if (!is.null(target) && (!is.numeric(target) || length(target) != 1 || !is.finite(target) ||
                           !is.null(names(target)))) {
    stop("For a numeric covariate, `target` is the target mean: a single unnamed number, as 45.")
  }
  if (is.null(resolution)) {
    resolution <- .rkn_default_resolution(x[known])
  } else if (!is.numeric(resolution) || length(resolution) != 1 || !is.finite(resolution) || resolution <= 0) {
    stop("`resolution` must be a single positive number.")
  }
  lowest <- min(x[known])
  steps <- round((x - lowest) / resolution)
  steps[!known] <- 0
  error <- abs(x[known] - lowest - resolution * steps[known])
  list(type = "mean", p = if (is.null(target)) mean(x[known]) else target, cohort = is.null(target),
       cat = steps, miss = as.numeric(!known), lowest = lowest, resolution = resolution,
       rounded = any(error > 1e-9 * pmax(1, abs(x[known]))))
}

#' Internal helper: default resolution of a numeric covariate
#'
#' Integer values spanning at most 100 keep a resolution of 1; otherwise, the smallest
#' of 1, 2 or 5 times a power of 10 that divides the range into at most 100 steps.
#'
#' @param v The known values.
#' @return The resolution.
#'
#' @keywords internal
#' @noRd
.rkn_default_resolution <- function(v) {
  span <- max(v) - min(v)
  if (span == 0 || (all(v == round(v)) && span <= 100)) return(1)
  for (power in floor(log10(span / 100)) + 0:1) {
    for (m in c(1, 2, 5)) {
      if (span / (m * 10^power) <= 100 + 1e-9) return(m * 10^power)
    }
  }
}

#' Internal helper: distance of the totals of a selection from the target
#'
#' @param tgt The result of \code{.rkn_global_target()}.
#' @return A function of the totals \code{C} (in the category, or steps of the values)
#'   and \code{M} (missing values) and of the size \code{S}, giving the distance of the
#'   proportion or mean among the individuals with a known value from the target, or
#'   \code{Inf} when none has one.
#'
#' @keywords internal
#' @noRd
.rkn_global_distance <- function(tgt) {
  p <- tgt$p
  lowest <- tgt$lowest
  resolution <- tgt$resolution
  if (tgt$type == "proportion") {
    function(C, M, S) {
      known <- S - M
      d <- abs(C - p * known) / known
      d[known <= 0] <- Inf
      d
    }
  } else {
    function(C, M, S) {
      known <- S - M
      d <- abs(lowest + resolution * C / known - p)
      d[known <= 0] <- Inf
      d
    }
  }
}
