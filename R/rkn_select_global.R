#' Selection with a target composition of the whole selection
#'
#' @description
#' Selects a maximum set of unrelated individuals whose composition, over the whole
#' selection, is as close as possible to a target: the proportion of one category of a
#' covariate, such as the proportion of women. Unlike \code{\link{rkn_select_criteria}},
#' which applies its criteria within each family, the families are chosen jointly: the
#' choice in one family depends on what the others contribute. For instance, with ten
#' unrelated women and ten related pairs of a woman and a man, choosing within each pair
#' gives fifteen women on average, while choosing the man in every pair gives ten women
#' and ten men.
#'
#' @inheritParams rkn_select_random
#' @param covariate Character. Name of a categorical vertex attribute of the network
#'   (character, factor or logical), such as \code{"sex"} for \code{V(g)$sex}.
#' @param target The target proportion of one category among the selected
#'   individuals, named by the category, as \code{c(F = 0.5)}; for a covariate with two
#'   categories, the proportions of both may be given, summing to 1, and the first is
#'   used. If \code{NULL} (default), the target is the proportion of the first category,
#'   in alphabetical order, in the whole network, so that the selection keeps the
#'   composition of the cohort; this requires a covariate with two categories.
#' @param max_loss Non-negative integer. Largest number of individuals that may be
#'   selected below the maximum, to bring the proportion closer to the target. With
#'   \code{0} (default) the selection keeps the maximum size. \code{\link{rkn_tradeoff}}
#'   shows, before selecting, the proportion that each number of individuals given up
#'   allows.
#'
#' @details
#' The selection is chosen in three steps: first the proportion closest to the target;
#' then, among the selections that reach it, the largest; then one of the selections
#' with that proportion and that size, drawn uniformly at random. With
#' \code{max_loss > 0}, individuals are given up even for a small improvement of the
#' proportion. The proportion is computed among the selected individuals with a known
#' value of the covariate; the number of selected individuals with a missing value is
#' reported by \code{print()}.
#'
#' The optimisation is exact. For each family, the exact algorithm described in the
#' section Search over maximum sets of \code{\link{rkn_selection}} counts its sets of
#' unrelated members by size and by number of members in the category; dynamic
#' programming then combines the families. The counts are kept on the logarithmic scale,
#' so that rare combinations, as for an extreme target, are not lost, and
#' \code{print()} reports the number of selections reaching the best proportion.
#'
#' The computation of a family is not attempted when it would exceed a budget of memory
#' and time, as for the largest families of genetic isolates at the third degree of
#' relationship. Such a family contributes the maximum sets with the smallest and the
#' largest number of members in the category, found exactly, and those of 100 uniform
#' draws; a family too intricate even for this is selected as by
#' \code{\link{rkn_select_random}}. The selection then keeps the maximum size, but the
#' proportion may not be the best possible and the draw is not uniform within these
#' families, whose number \code{print()} reports. If the tables that combine the
#' families would exceed about 80 MB, which requires many individuals with a missing
#' value or a large \code{max_loss}, the selection stops with an error before
#' allocating them.
#'
#' @return An object of class \code{\link{rkn_selection}} with \code{method = "global"};
#'   its element \code{global} describes the composition reached (see
#'   \code{\link{rkn_selection}}).
#'
#' @seealso \code{\link{rkn_tradeoff}} for the proportions reachable by each number of
#'   individuals given up; \code{\link{rkn_selection}}, \code{\link{rkn_select_criteria}},
#'   \code{\link{rkn_select_random}}.
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
#' @export
rkn_select_global <- function(g, covariate, target = NULL, max_loss = 0, seed = NULL, diagnostics = TRUE) {

  if (missing(covariate)) covariate <- NULL
  .rkn_global_check_input(g, covariate)
  if (!is.numeric(max_loss) || length(max_loss) != 1 || is.na(max_loss) || max_loss < 0 ||
      max_loss != round(max_loss)) {
    stop("`max_loss` must be a single non-negative integer.")
  }
  tgt <- .rkn_global_target(igraph::vertex_attr(g, covariate), target)

  if (!is.null(seed)) {
    withr::local_preserve_seed()
    set.seed(seed)
  }

  ids <- igraph::V(g)$name
  fam <- .rkn_global_families(g, tgt, max_loss)
  comp <- fam$comp

  p <- tgt$p
  distance <- function(C, M, S) {
    known <- S - M
    d <- abs(C - p * known) / known
    d[known <= 0] <- Inf
    d
  }
  # The network cannot give up more individuals than its maximum size.
  res <- .rkn_global_select(fam$fams, min(max_loss, fam$max_size), distance)
  sets <- .rkn_global_realise(fam$fams, res$choice)
  selected <- rep(FALSE, length(ids))
  for (i in seq_along(sets)) selected[fam$members[[i]][sets[[i]]]] <- TRUE

  global <- list(
    covariate = covariate, category = tgt$category, target = p, cohort = tgt$cohort,
    proportion = res$c / (res$size - res$m), n_category = res$c, n_known = res$size - res$m, n_missing = res$m,
    size = res$size, max_size = res$size + res$j, max_loss = max_loss, log_count = res$log_count,
    exact = res$exact, n_inexact = sum(!fam$exact)
  )

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
      n_discarded = tabulate(family_id[!selected], nbins = comp$no),
      n_category = tabulate(family_id[selected & tgt$cat == 1], nbins = comp$no),
      n_missing = tabulate(family_id[selected & tgt$miss == 1], nbins = comp$no)
    )
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
#' Groups the individuals without relatives by kind and computes the options of every
#' family (\code{.rkn_global_options()}). A family cannot give up more members than it
#' has, so that its options are computed for at most that many: a larger \code{k} would
#' only enlarge its table.
#'
#' @param g The network.
#' @param tgt The result of \code{.rkn_global_target()}.
#' @param k Largest number of individuals below the maximum size of the network.
#' @param previous A result of this function for another \code{k}: the options of the
#'   families for which the limit is unchanged are reused.
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
    if (!is.null(previous) && previous$cap[f] == cap[f]) {
      options[[f]] <- previous$options[[f]]
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

#' Internal helper: target category and proportion of a global selection
#'
#' @param x The covariate, one value per individual.
#' @param target The \code{target} argument of \code{rkn_select_global()}.
#' @return A list with \code{category}, \code{p}, \code{cohort} (whether the target is
#'   the proportion of the network), and the 0/1 vectors \code{cat} (in the category)
#'   and \code{miss} (missing value).
#'
#' @keywords internal
#' @noRd
.rkn_global_target <- function(x, target) {
  if (is.numeric(x)) {
    stop("rkn_select_global() supports categorical covariates (character, factor or logical); ",
         "store the categories as character or factor.")
  }
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
  list(category = category, p = p, cohort = is.null(target),
       cat = as.numeric(known & x == category), miss = as.numeric(!known))
}
