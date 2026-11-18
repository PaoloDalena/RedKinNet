runif_after <- function(f) {
  set.seed(123)
  f()
  stats::runif(3)
}

reference <- function() {
  set.seed(123)
  stats::runif(3)
}

test_that("seeded network-level selections leave the caller's random number stream untouched", {
  g <- make_network(41)
  igraph::V(g)$sex <- rep(c("F", "M"), length.out = igraph::vcount(g))
  crit <- rkn_criterion_level("sex", "F")
  expect_identical(runif_after(function() rkn_select_random(g, seed = 1)), reference())
  expect_identical(runif_after(function() rkn_select_criteria(g, crit, seed = 1)), reference())
})

test_that("seeded per-family selections leave the caller's random number stream untouched", {
  g5 <- igraph::make_ring(5)
  igraph::V(g5)$name <- letters[1:5]
  g7 <- igraph::make_ring(7)
  igraph::V(g7)$name <- letters[1:7]
  fd5 <- fam_df_from_graph(g5)
  fd7 <- fam_df_from_graph(g7)
  fd5$sex <- c("F", "M", "F", "M", "F")
  fd7$sex <- c("F", "M", "F", "M", "F", "M", "F")
  crit <- rkn_criterion_level("sex", "F")
  expect_identical(runif_after(function() .rkn_select_random_upto5(fd5, seed = 1, g_sub = g5)), reference())
  expect_identical(runif_after(function() .rkn_select_random_6plus(fd7, g7, seed = 1)), reference())
  expect_identical(runif_after(function() .rkn_select_criteria_upto5(fd5, crit, seed = 1, g_sub = g5)), reference())
  expect_identical(runif_after(function() .rkn_select_criteria_6plus(fd7, g7, crit, seed = 1)), reference())
})

test_that("a session without a random seed stays without one after a seeded selection", {
  g <- make_network(42)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  rkn_select_random(g, seed = 1)
  expect_false(exists(".Random.seed", envir = globalenv()))
})

test_that("the plot layout leaves the caller's random number stream untouched", {
  g <- make_network(41)
  expect_identical(runif_after(function() .rkn_layout(g, "components", seed = 1)), reference())
  expect_identical(runif_after(function() .rkn_layout(g, "grid", seed = 1)), reference())
})

test_that("seeded uniform selections leave the caller's random number stream untouched", {
  g <- make_network(43)
  expect_identical(runif_after(function() rkn_select_random(g, seed = 1, uniform = TRUE)), reference())
})

test_that("consecutive seeds give independent selections", {
  # 300 identical pairs: if family i used the seed seed + i, the selection with seed 2
  # would repeat in pair i the choice made with seed 1 in pair i + 1.
  n <- 300
  g <- igraph::make_empty_graph(2 * n, directed = FALSE)
  g <- igraph::add_edges(g, as.vector(rbind(seq(1, 2 * n, 2), seq(2, 2 * n, 2))))
  igraph::V(g)$name <- paste0("v", seq_len(2 * n))
  igraph::V(g)$sex <- "F"
  first_kept <- function(sel) seq(1, 2 * n, 2) %in% which(igraph::V(g)$name %in% sel$kept)
  selectors <- list(
    function(seed) rkn_select_random(g, seed = seed, diagnostics = FALSE),
    function(seed) rkn_select_random(g, seed = seed, uniform = TRUE, diagnostics = FALSE),
    function(seed) rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = seed, diagnostics = FALSE)
  )
  for (select in selectors) {
    a <- first_kept(select(1))
    b <- first_kept(select(2))
    expect_lt(mean(b[-n] == a[-1]), 0.75)
    expect_lt(mean(b == a), 0.75)
  }
})
