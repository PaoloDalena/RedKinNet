adjacency <- function(g) as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE)) > 0

draw_sets <- function(sol, ids, n_draws) {
  drawn <- character(n_draws)
  for (d in seq_len(n_draws)) drawn[d] <- paste(sort(ids[.rkn_mis_sample(sol)]), collapse = ",")
  drawn
}

test_that("the exact engine counts the maximum independent sets of every connected graph up to 7 nodes", {
  pool <- thin_on_cran(atlas_connected(1, 7))
  problems <- character(0)
  for (nm in names(pool)) {
    g <- pool[[nm]]
    ref <- brute_force_mis(g)
    sol <- .rkn_mis_exact(adjacency(g))
    if (sol$alpha != ref$alpha || round(exp(sol$log_count)) != length(ref$sets)) {
      problems <- c(problems, nm)
    }
  }
  expect_identical(problems, character(0))
})

test_that("the exact engine agrees with brute force on random and disconnected graphs", {
  set.seed(5)
  for (k in 1:30) {
    n <- sample(8:14, 1)
    g <- igraph::sample_gnp(n, runif(1, 0.1, 0.5))
    igraph::V(g)$name <- as.character(seq_len(n))
    ref <- brute_force_mis(g)
    sol <- .rkn_mis_exact(adjacency(g))
    expect_equal(sol$alpha, ref$alpha)
    expect_equal(round(exp(sol$log_count)), length(ref$sets))
  }
})

test_that("draws are maximum independent sets", {
  set.seed(6)
  for (k in 1:10) {
    g <- igraph::sample_gnp(40, 0.08)
    A <- adjacency(g)
    sol <- .rkn_mis_exact(A, max_width = Inf)
    for (d in 1:5) {
      x <- .rkn_mis_sample(sol)
      expect_identical(sum(x), as.integer(sol$alpha))
      expect_false(any(A[x, x]))
    }
  }
})

test_that("draws are uniform over the maximum independent sets", {
  set.seed(7)
  path6 <- igraph::make_ring(6, circular = FALSE)
  drawn <- draw_sets(.rkn_mis_exact(adjacency(path6)), as.character(1:6), 4000)
  expect_setequal(unique(drawn), c("1,3,5", "1,3,6", "1,4,6", "2,4,6"))
  expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)

  g <- igraph::sample_gnp(12, 0.3)
  ids <- as.character(seq_len(12))
  igraph::V(g)$name <- ids
  ref <- brute_force_mis(g)
  expect_gt(length(ref$sets), 2)
  drawn <- draw_sets(.rkn_mis_exact(adjacency(g)), ids, 200 * length(ref$sets))
  expect_setequal(unique(drawn), vapply(ref$sets, function(s) paste(sort(s), collapse = ","), ""))
  expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
})

test_that("the width of the elimination ordering reflects the structure of the graph", {
  expect_identical(.rkn_elimination_order(adjacency(igraph::make_tree(30, 2, mode = "undirected")))$width, 1L)
  expect_identical(.rkn_elimination_order(adjacency(igraph::make_full_graph(8)))$width, 7L)
  expect_identical(.rkn_elimination_order(adjacency(igraph::make_ring(10)))$width, 2L)
  sol <- .rkn_mis_exact(adjacency(igraph::make_full_graph(8)), max_width = 5)
  expect_false(sol$feasible)
  expect_identical(sol$width, 7L)
})

test_that("very large numbers of maximum independent sets are kept on the log scale", {
  pairs <- igraph::make_graph(seq_len(2 * 2000), directed = FALSE)
  sol <- .rkn_mis_exact(adjacency(pairs))
  expect_identical(sol$alpha, 2000)
  expect_equal(sol$log_count, 2000 * log(2))
})

test_that("dominated vertices are removed without reducing the maximum size", {
  set.seed(8)
  for (k in 1:20) {
    g <- igraph::sample_gnp(sample(8:14, 1), runif(1, 0.15, 0.5))
    A <- adjacency(g)
    keep <- .rkn_dominance_reduce(A)
    expect_equal(.rkn_mis_exact(A[keep, keep, drop = FALSE])$alpha, .rkn_mis_exact(A)$alpha)
  }
  # A hub related to everyone in a chain is dominated by the end of the chain.
  fan <- igraph::make_ring(8, circular = FALSE) + igraph::vertex(1)
  fan <- igraph::add_edges(fan, as.vector(rbind(9, 1:8)))
  expect_false(.rkn_dominance_reduce(adjacency(fan))[9])
})

test_that("draws fall back on the dominance reduction and then stop when a graph is too intricate", {
  fan <- igraph::make_ring(8, circular = FALSE) + igraph::vertex(1)
  fan <- igraph::add_edges(fan, as.vector(rbind(9, 1:8)))
  igraph::V(fan)$name <- as.character(1:9)
  drawn <- .rkn_draw_maximum_set(fan, max_width = 1)
  expect_identical(drawn$method, "exact_reduced")
  expect_length(drawn$kept, 4)
  expect_false(any(adjacency(fan)[drawn$kept, drawn$kept]))

  cycle <- igraph::make_ring(7)
  igraph::V(cycle)$name <- as.character(1:7)
  expect_identical(.rkn_draw_maximum_set(cycle)$method, "exact")
  expect_error(.rkn_draw_maximum_set(cycle, max_width = 1), "too intricate for an exact selection")
})

test_that("all maximum independent sets are listed once", {
  set.seed(9)
  for (k in 1:20) {
    g <- igraph::sample_gnp(sample(6:13, 1), runif(1, 0.1, 0.5))
    ids <- as.character(seq_len(igraph::vcount(g)))
    igraph::V(g)$name <- ids
    ref <- brute_force_mis(g)
    sets <- .rkn_mis_enumerate(.rkn_mis_exact(adjacency(g)), limit = 1e4)
    listed <- character(ncol(sets))
    for (j in seq_len(ncol(sets))) listed[j] <- paste(sort(ids[sets[, j]]), collapse = ",")
    expect_false(anyDuplicated(listed) > 0)
    expect_setequal(listed, vapply(ref$sets, function(s) paste(sort(s), collapse = ","), ""))
  }
  pairs <- igraph::make_graph(seq_len(2 * 20), directed = FALSE)
  expect_null(.rkn_mis_enumerate(.rkn_mis_exact(adjacency(pairs)), limit = 1e4))
})

test_that("many draws at once are uniform maximum independent sets", {
  set.seed(10)
  path6 <- igraph::make_ring(6, circular = FALSE)
  x <- .rkn_mis_sample_many(.rkn_mis_exact(adjacency(path6)), 4000)
  drawn <- character(ncol(x))
  for (j in seq_len(ncol(x))) drawn[j] <- paste(which(x[, j]), collapse = ",")
  expect_setequal(unique(drawn), c("1,3,5", "1,3,6", "1,4,6", "2,4,6"))
  expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
})

test_that("weights rank the maximum sets by an indicator after their size", {
  path6 <- igraph::make_ring(6, circular = FALSE)
  # Preferring vertices 2 and 5 among the maximum sets of size 3 leaves {2,4,6} and {1,3,5}
  # with one preferred vertex each, and {1,3,6}, {1,4,6} with none.
  hit <- c(0, 1, 0, 0, 1, 0)
  sol <- .rkn_mis_exact(adjacency(path6), weight = 7 + hit)
  expect_equal(sol$alpha, 3 * 7 + 1)
  sets <- .rkn_mis_enumerate(sol, limit = 100)
  listed <- character(ncol(sets))
  for (j in seq_len(ncol(sets))) listed[j] <- paste(which(sets[, j]), collapse = ",")
  expect_setequal(listed, c("1,3,5", "2,4,6"))
})

test_that("inclusion probabilities match the share of maximum independent sets containing each vertex", {
  check <- function(g) {
    ids <- as.character(seq_len(igraph::vcount(g)))
    igraph::V(g)$name <- ids
    sets <- brute_force_mis(g)$sets
    expected <- numeric(length(ids))
    for (s in sets) expected <- expected + (ids %in% s)
    all.equal(.rkn_mis_inclusion(.rkn_mis_exact(adjacency(g))), expected / length(sets))
  }
  pool <- thin_on_cran(atlas_connected(1, 7))
  problems <- character(0)
  for (nm in names(pool)) if (!isTRUE(check(pool[[nm]]))) problems <- c(problems, nm)
  expect_identical(problems, character(0))
  set.seed(12)
  for (k in 1:30) {
    g <- igraph::sample_gnp(sample(8:14, 1), runif(1, 0.1, 0.5))
    expect_true(isTRUE(check(g)))
  }
  path6 <- .rkn_mis_inclusion(.rkn_mis_exact(adjacency(igraph::make_ring(6, circular = FALSE))))
  expect_equal(path6, c(0.75, 0.25, 0.5, 0.5, 0.25, 0.75))
})
