adjacency <- function(g) as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE)) > 0

proportion_distance <- function(p) {
  function(C, M, S) {
    D <- S - M
    d <- abs(C - p * D) / D
    d[D <= 0] <- Inf
    d
  }
}

# Global selection of a network by proportion of the category, with the single
# individuals grouped as in the package.
global_selection <- function(g, cat, miss, k, p, ...) {
  comp <- igraph::components(g)
  single <- which(comp$csize[comp$membership] == 1)
  fams <- list()
  idx <- list()
  for (grp in .rkn_global_singletons(cat[single], miss[single], k)) {
    fams[[length(fams) + 1]] <- grp
    idx[[length(idx) + 1]] <- single[grp$members]
  }
  for (f in which(comp$csize > 1)) {
    v <- which(comp$membership == f)
    fams[[length(fams) + 1]] <- .rkn_global_options(adjacency(igraph::induced_subgraph(g, v)), cat[v], miss[v], k, ...)
    idx[[length(idx) + 1]] <- v
  }
  res <- .rkn_global_select(fams, k, proportion_distance(p))
  sets <- .rkn_global_realise(fams, res$choice)
  sel <- rep(FALSE, igraph::vcount(g))
  for (f in seq_along(sets)) sel[idx[[f]][sets[[f]]]] <- TRUE
  list(res = res, sel = sel)
}

# All the selections with the best proportion and, among them, the largest size.
best_selections <- function(g, cat, miss, k, p) {
  A <- adjacency(g)
  n <- nrow(A)
  labels <- character(0)
  d <- numeric(0)
  size <- integer(0)
  for (s in 0:(2^n - 1)) {
    sel <- bitwAnd(s, 2^(seq_len(n) - 1)) > 0
    if (sum(sel) > 1 && any(A[sel, sel])) next
    labels <- c(labels, paste(which(sel), collapse = ","))
    size <- c(size, sum(sel))
    D <- sum(sel) - sum(miss[sel])
    d <- c(d, if (D > 0) abs(sum(cat[sel]) - p * D) / D else Inf)
  }
  keep <- size >= max(size) - k
  labels <- labels[keep]
  d <- d[keep]
  size <- size[keep]
  best <- d <= min(d) + 1e-12
  labels[best & size == max(size[best])]
}

random_network <- function(n_families) {
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1), atlas_connected(1, 1))
  repeat {
    parts <- pool[sample(length(pool), n_families)]
    for (i in seq_along(parts)) igraph::V(parts[[i]])$name <- paste0(i, "_", seq_len(igraph::vcount(parts[[i]])))
    g <- do.call(igraph::disjoint_union, unname(parts))
    if (igraph::vcount(g) <= 12) return(g)
  }
}

test_that("global selections reach the best proportion with the largest size, and count the selections that do", {
  set.seed(19)
  for (i in 1:30) {
    g <- random_network(sample(2:4, 1))
    n <- igraph::vcount(g)
    miss <- stats::rbinom(n, 1, 0.15)
    cat <- sample(0:1, n, replace = TRUE) * (1 - miss)
    k <- sample(0:2, 1)
    p <- stats::runif(1)
    best <- best_selections(g, cat, miss, k, p)
    out <- global_selection(g, cat, miss, k, p)
    expect_true(paste(which(out$sel), collapse = ",") %in% best)
    expect_equal(exp(out$res$log_count), length(best))
    expect_true(out$res$exact)
  }
})

test_that("global selections are drawn uniformly among the best selections", {
  set.seed(20)
  tested <- 0
  while (tested < 3) {
    g <- random_network(3)
    n <- igraph::vcount(g)
    miss <- stats::rbinom(n, 1, 0.15)
    cat <- sample(0:1, n, replace = TRUE) * (1 - miss)
    k <- sample(0:1, 1)
    p <- stats::runif(1)
    best <- best_selections(g, cat, miss, k, p)
    if (length(best) < 2 || length(best) > 8) next
    tested <- tested + 1
    drawn <- character(150 * length(best))
    for (d in seq_along(drawn)) drawn[d] <- paste(which(global_selection(g, cat, miss, k, p)$sel), collapse = ",")
    expect_setequal(unique(drawn), best)
    expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
  }
})

test_that("extreme targets are reached among more selections than double precision can count", {
  # 1100 related pairs of a woman and a man: 2^1100 maximum sets, one with all women,
  # choose(1100, 550) with half of them.
  pairs <- igraph::make_graph(seq_len(2200), directed = FALSE)
  cat <- rep(c(1, 0), 1100)
  all_women <- global_selection(pairs, cat, rep(0, 2200), 0, 1)
  expect_identical(sum(cat[all_women$sel]), 1100)
  expect_equal(all_women$res$log_count, 0)
  half <- global_selection(pairs, cat, rep(0, 2200), 0, 0.5)
  expect_identical(sum(cat[half$sel]), 550)
  expect_equal(half$res$log_count, lchoose(1100, 550))
})

test_that("families beyond the budget or too intricate are included with the options found, and flagged", {
  set.seed(21)
  petersen <- igraph::make_graph("Petersen")
  cat <- rep(c(1, 0), 5)
  partial <- global_selection(petersen, cat, rep(0, 10), 0, 0.5, max_cells = 1)
  expect_false(partial$res$exact)
  expect_true(is.na(partial$res$log_count))
  expect_identical(sum(partial$sel), 4L)
  expect_identical(igraph::ecount(igraph::induced_subgraph(petersen, which(partial$sel))), 0)

  ring <- igraph::make_ring(7)
  fallback <- function() c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)
  fixed <- .rkn_global_options(adjacency(ring), rep(1, 7), rep(0, 7), 0, fallback = fallback, max_cells = 1, max_width = 1)
  expect_false(fixed$exact)
  expect_identical(fixed$options$c, 3)
  expect_error(.rkn_global_options(adjacency(ring), rep(1, 7), rep(0, 7), 0, max_cells = 1, max_width = 1),
               "too intricate")
})

test_that("the memory of the global selection is checked before the tables are built", {
  b <- .rkn_global_singletons(rep(c(1, 0), 50), rep(0, 100), 2)
  expect_error(.rkn_global_select(b, 2, proportion_distance(0.5), max_cells = 10), "table cells")
})
