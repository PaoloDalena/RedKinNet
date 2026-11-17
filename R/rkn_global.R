#' Internal helper: options of a family for a global selection
#'
#' Lists what a family can contribute to a global selection: for each number \code{j}
#' of individuals below its maximum size (0 to \code{k}), each combination of
#' \code{c}, the number of selected members in the category, and \code{m}, the number of
#' selected members with a missing value, with \code{lc}, the logarithm of the number of
#' sets of unrelated members that give it. The two counts are combined into one value,
#' \code{c + (n_c + 1) m} with \code{n_c} the members in the category, so that the profile
#' of \code{.rkn_family_profile()} counts both at once. Single individuals and pairs are
#' listed directly. A family beyond the budget of the profile contributes the maximum
#' sets of its partial profile, a family too intricate even for that the set returned
#' by \code{fallback}; in both cases the options are not exhaustive, their counts are
#' not known, and the family is marked \code{exact = FALSE}.
#'
#' @param A Symmetric logical adjacency matrix of a connected family.
#' @param cat,miss 0/1 vectors: member in the category, member with a missing value.
#' @param k Largest number of individuals below the maximum size.
#' @param fallback A function returning a maximum set of the family as a logical vector,
#'   or \code{NULL}.
#' @param max_cells,max_work Budget of the profile, as in \code{.rkn_family_profile()}.
#' @param max_width Largest elimination width of the partial profile.
#' @return A list with \code{alpha}, \code{exact}, \code{options} (a list of the vectors
#'   \code{j}, \code{c}, \code{m}, \code{lc} and \code{v}, the combined value), and
#'   either \code{A}, \code{x}, \code{elim} and \code{k} (options drawn from the
#'   profile) or \code{sets} (for each option, the list of the sets that give it).
#'
#' @keywords internal
#' @noRd
.rkn_global_options <- function(A, cat, miss, k, fallback = NULL, max_cells = 4e6, max_work = 4e9,
                                max_width = 20) {
  n <- nrow(A)
  listed <- function(sets, alpha, exact) {
    j <- c_ <- m_ <- numeric(length(sets))
    key <- character(length(sets))
    for (i in seq_along(sets)) {
      j[i] <- alpha - sum(sets[[i]])
      c_[i] <- sum(cat[sets[[i]]])
      m_[i] <- sum(miss[sets[[i]]])
      key[i] <- paste(j[i], c_[i], m_[i])
    }
    keys <- unique(key)
    first <- integer(length(keys))
    size <- integer(length(keys))
    option_sets <- vector("list", length(keys))
    for (g in seq_along(keys)) {
      members <- which(key == keys[g])
      first[g] <- members[1]
      size[g] <- length(members)
      option_sets[[g]] <- sets[members]
    }
    list(alpha = alpha, exact = exact,
         options = list(j = j[first], c = c_[first], m = m_[first],
                        lc = if (exact) log(size) else rep(0, length(keys)), v = rep(NA_real_, length(keys))),
         sets = option_sets)
  }

  if (n <= 2) {
    sets <- if (n == 1) list(TRUE) else list(c(TRUE, FALSE), c(FALSE, TRUE))
    if (k >= 1) sets <- c(sets, list(rep(FALSE, n)))
    return(listed(sets, 1, TRUE))
  }

  base <- sum(cat) + 1
  x <- cat + base * miss
  elim <- .rkn_elimination_order(A)
  prof <- .rkn_family_profile(A, x, k, elim, max_cells, max_work)
  if (!is.null(prof)) {
    nz <- which(prof$counts > 0, arr.ind = TRUE)
    v <- nz[, 2] - 1
    return(list(alpha = prof$alpha, exact = TRUE,
                options = list(j = nz[, 1] - 1, c = v %% base, m = v %/% base,
                               lc = log(prof$counts[nz]) + prof$log_scale, v = v),
                A = A, x = x, elim = elim, k = k))
  }
  part <- .rkn_partial_profile(A, cat, max_width = max_width)
  if (!is.null(part)) return(listed(part$sets, part$alpha, FALSE))
  if (is.null(fallback)) stop("A family is too intricate for a global selection.", call. = FALSE)
  set <- fallback()
  listed(list(set), sum(set), FALSE)
}

#' Internal helper: single individuals of a global selection, grouped
#'
#' Single individuals of the same kind (in the category, not in it, with a missing
#' value) are interchangeable, so each kind forms one group whose options are the
#' number \code{t} of members left out, from 0 to \code{k}, in \code{choose(n, t)}
#' ways. This keeps the number of families to combine small when \code{k > 0}, where
#' every single individual could otherwise be left out on its own.
#'
#' @param cat,miss 0/1 vectors over the single individuals.
#' @param k Largest number of individuals below the maximum size.
#' @return A list of groups, each in the format of \code{.rkn_global_options()} with
#'   \code{group}, the number of members, and \code{members}, their positions in
#'   \code{cat}.
#'
#' @keywords internal
#' @noRd
.rkn_global_singletons <- function(cat, miss, k) {
  groups <- list()
  for (kind in list(c(1, 0), c(0, 0), c(0, 1))) {
    members <- which(cat == kind[1] & miss == kind[2])
    n <- length(members)
    if (n == 0) next
    t <- 0:min(k, n)
    groups[[length(groups) + 1]] <- list(
      alpha = n, exact = TRUE, group = n, members = members,
      options = list(j = t, c = (n - t) * kind[1], m = (n - t) * kind[2], lc = lchoose(n, t),
                     v = rep(NA_real_, length(t)))
    )
  }
  groups
}

#' Internal helper: logarithm of the sum of two exponentials, elementwise
#'
#' @keywords internal
#' @noRd
.rkn_log_add <- function(a, b) {
  hi <- pmax(a, b)
  out <- hi + log1p(exp(-abs(a - b)))
  out[hi == -Inf] <- -Inf
  out
}

#' Internal helper: best global selection over the options of the families
#'
#' Combines the options of the families (\code{.rkn_global_options()}) into a table
#' over the totals of the selection: the number \code{J} of individuals below the
#' maximum size of the network, the number \code{C} of selected individuals in the
#' category and the number \code{M} with a missing value. Each entry holds the logarithm
#' of the number of selections reaching it, so that combinations that are rare but
#' optimal, as for an extreme target, are not lost to rounding. The best totals are
#' those at the smallest \code{distance}; among them, those with the largest selection;
#' one of the selections reaching them is drawn uniformly at random, by going back
#' through the families and choosing the option of each in proportion to the number of
#' selections of the families before it that complete it. Each table covers only the
#' totals reachable by the families combined so far, and only the tables at regular
#' intervals are kept, the others being recomputed when going back, so that the memory
#' grows with the square root of the number of families. The memory needed is computed
#' before anything is allocated.
#'
#' @param fams A list of results of \code{.rkn_global_options()}.
#' @param k Largest number of individuals below the maximum size.
#' @param distance A function of the totals \code{C}, \code{M} and of the size \code{S},
#'   arrays of the same shape, giving the distance from the target (\code{Inf} where it
#'   is not defined).
#' @param max_cells Largest number of cells of the tables kept at once. Beyond it,
#'   the function stops with an error of class \code{rkn_table_limit}.
#' @param draw If \code{FALSE}, only the table of the totals is computed, without
#'   keeping the tables needed to go back, and \code{distance} is not used.
#' @return A list with \code{choice}, the option chosen in each family; the totals
#'   \code{size}, \code{j}, \code{c} and \code{m} of the selection; \code{distance};
#'   \code{log_count}, the logarithm of the number of selections reaching the best
#'   totals (\code{NA} when some family is not \code{exact}); \code{exact}; and
#'   \code{table}, the table of the totals. With \code{draw = FALSE}, a list with
#'   \code{table}, \code{alpha} (the maximum size), \code{c_base} and \code{m_base}
#'   (the totals of the first cell of the table) and \code{exact}.
#'
#' @keywords internal
#' @noRd
.rkn_global_select <- function(fams, k, distance, max_cells = 1e7, draw = TRUE) {
  n_fam <- length(fams)
  n_opt <- integer(n_fam)
  alpha <- 0
  c_base <- 0
  m_base <- 0
  c_width <- numeric(n_fam)
  m_width <- numeric(n_fam)
  exact <- TRUE
  # The tables hold the totals above the sum of the smallest counts of the families.
  for (f in seq_len(n_fam)) {
    o <- fams[[f]]$options
    n_opt[f] <- length(o$j)
    alpha <- alpha + fams[[f]]$alpha
    c_base <- c_base + min(o$c)
    m_base <- m_base + min(o$m)
    c_width[f] <- max(o$c) - min(o$c)
    m_width[f] <- max(o$m) - min(o$m)
    fams[[f]]$options$c <- o$c - min(o$c)
    fams[[f]]$options$m <- o$m - min(o$m)
    exact <- exact && fams[[f]]$exact
  }
  # The table grows with the families combined so far: families with missing values,
  # which widen the third dimension, come last.
  free <- which(n_opt > 1)
  free <- free[order(m_width[free] > 0, c_width[free], m_width[free])]
  n_free <- length(free)
  block <- max(1, ceiling(sqrt(n_free)))
  box <- (k + 1) * (cumsum(c(0, c_width[free])) + 1) * (cumsum(c(0, m_width[free])) + 1)
  starts <- if (n_free > 0) seq(1, n_free, by = block) else integer(0)
  kept <- box[c(starts, n_free + 1)]
  peak <- 0
  for (first in starts) peak <- max(peak, sum(box[first:min(first + block - 1, n_free)]))
  # Without going back, only two consecutive tables are held at once.
  need <- if (draw) sum(kept) + peak else max(box[-1] + box[-length(box)], box[1])
  if (need > max_cells) {
    msg <- sprintf(paste0("The global selection would keep %.1e table cells (limit %.1e). Reduce max_loss, or ",
                          "the number of individuals with a missing value of the covariate."), need, max_cells)
    stop(structure(class = c("rkn_table_limit", "error", "condition"), list(message = msg, call = NULL)))
  }

  combine <- function(tab, o) {
    d <- dim(tab)
    out <- array(-Inf, c(k + 1, d[2] + max(o$c), d[3] + max(o$m)))
    for (i in seq_along(o$j)) {
      if (o$j[i] > k) next
      sj <- seq_len(k + 1 - o$j[i])
      tj <- o$j[i] + sj
      tc <- o$c[i] + seq_len(d[2])
      tm <- o$m[i] + seq_len(d[3])
      out[tj, tc, tm] <- .rkn_log_add(out[tj, tc, tm, drop = FALSE], tab[sj, , , drop = FALSE] + o$lc[i])
    }
    out
  }

  # Families with a single option shift the starting point.
  start <- 0
  log_start <- 0
  for (f in which(n_opt == 1)) {
    start <- start + fams[[f]]$options$j
    log_start <- log_start + fams[[f]]$options$lc
  }
  tab <- array(-Inf, c(k + 1, 1, 1))
  tab[start + 1, 1, 1] <- log_start
  checkpoints <- list(tab)
  for (i in seq_len(n_free)) {
    tab <- combine(tab, fams[[free[i]]]$options)
    if (draw && i %% block == 0 && i < n_free) checkpoints[[length(checkpoints) + 1]] <- tab
  }
  if (!draw) return(list(table = tab, alpha = alpha, c_base = c_base, m_base = m_base, exact = exact))

  size <- alpha - (slice.index(tab, 1) - 1)
  d <- distance(c_base + slice.index(tab, 2) - 1, m_base + slice.index(tab, 3) - 1, size)
  d[tab == -Inf] <- Inf
  d_min <- min(d)
  if (!is.finite(d_min)) stop("No selection has a defined value of the criterion.", call. = FALSE)
  best <- d <= d_min + 1e-12
  best <- best & size == max(size[best])
  lw <- tab[best]
  log_count <- max(lw) + log(sum(exp(lw - max(lw))))
  pick <- which(best)[sample.int(length(lw), 1, prob = exp(lw - max(lw)))]
  target <- as.vector(arrayInd(pick, dim(tab))) - 1
  totals <- target

  choice <- rep(1L, n_fam)
  for (b in rev(seq_along(checkpoints))) {
    first <- (b - 1) * block + 1
    last <- min(b * block, n_free)
    if (first > last) next
    before <- vector("list", last - first + 1)
    before[[1]] <- checkpoints[[b]]
    if (last > first) {
      for (i in first:(last - 1)) {
        before[[i - first + 2]] <- combine(before[[i - first + 1]], fams[[free[i]]]$options)
      }
    }
    for (i in last:first) {
      o <- fams[[free[i]]]$options
      prev <- before[[i - first + 1]]
      w <- rep(-Inf, length(o$j))
      for (q in seq_along(o$j)) {
        rest <- target - c(o$j[q], o$c[q], o$m[q])
        if (all(rest >= 0) && all(rest < dim(prev))) w[q] <- prev[rest[1] + 1, rest[2] + 1, rest[3] + 1] + o$lc[q]
      }
      q <- sample.int(length(w), 1, prob = exp(w - max(w)))
      choice[free[i]] <- q
      target <- target - c(o$j[q], o$c[q], o$m[q])
    }
  }
  stopifnot(all(target == c(start, 0, 0)))

  list(choice = choice, size = alpha - totals[1], j = totals[1], c = c_base + totals[2], m = m_base + totals[3],
       distance = d_min, log_count = if (exact) log_count else NA_real_, exact = exact, table = tab)
}

#' Internal helper: the sets of unrelated individuals chosen in each family
#'
#' Draws, in each family, one of the sets that give the option chosen by
#' \code{.rkn_global_select()}, uniformly at random: from the profile of the family,
#' recomputed with its messages, among the listed sets, or, for a group of single
#' individuals, by leaving out the chosen number of them.
#'
#' @param fams A list of results of \code{.rkn_global_options()}.
#' @param choice The option chosen in each family.
#' @return A list with one logical vector per family.
#'
#' @keywords internal
#' @noRd
.rkn_global_realise <- function(fams, choice) {
  out <- vector("list", length(fams))
  for (f in seq_along(fams)) {
    fam <- fams[[f]]
    q <- choice[f]
    if (!is.null(fam$group)) {
      out[[f]] <- rep(TRUE, fam$group)
      out[[f]][sample.int(fam$group, fam$options$j[q])] <- FALSE
    } else if (is.null(fam$sets)) {
      prof <- .rkn_family_profile(fam$A, fam$x, fam$k, fam$elim, max_cells = Inf, max_work = Inf, keep = TRUE)
      out[[f]] <- .rkn_profile_sample(prof, fam$options$j[q], fam$options$v[q])
    } else {
      sets <- fam$sets[[q]]
      out[[f]] <- sets[[sample.int(length(sets), 1)]]
    }
  }
  out
}
