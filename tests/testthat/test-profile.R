adjacency <- function(g) as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE)) > 0

# Counts of the independent sets of size alpha - j (rows, j = 0..k) by value sum
# (columns, from 0), by enumeration of all subsets.
profile_brute <- function(A, x, k) {
  n <- nrow(A)
  sizes <- integer(0)
  sums <- numeric(0)
  for (s in 0:(2^n - 1)) {
    sel <- which(bitwAnd(s, 2^(seq_len(n) - 1)) > 0)
    if (length(sel) > 1 && any(A[sel, sel])) next
    sizes <- c(sizes, length(sel))
    sums <- c(sums, sum(x[sel]))
  }
  alpha <- max(sizes)
  counts <- matrix(0, k + 1, sum(x) + 1)
  for (i in which(sizes >= alpha - k)) {
    counts[alpha - sizes[i] + 1, sums[i] + 1] <- counts[alpha - sizes[i] + 1, sums[i] + 1] + 1
  }
  list(alpha = alpha, counts = counts)
}

same_profile <- function(A, x, k) {
  p <- .rkn_family_profile(A, x, k)
  b <- profile_brute(A, x, k)
  p$alpha == b$alpha && isTRUE(all.equal(p$counts * exp(p$log_scale), b$counts))
}

# A hub related to both members of m pairs: 2^m maximum sets, one member of each pair.
windmill <- function(m) {
  edges <- integer(0)
  for (i in seq_len(m)) edges <- c(edges, 1, 2 * i, 1, 2 * i + 1, 2 * i, 2 * i + 1)
  igraph::make_graph(edges, directed = FALSE)
}

test_that("profiles count the independent sets of every connected graph up to 7 nodes by size and value", {
  pool <- thin_on_cran(atlas_connected(2, 7))
  set.seed(15)
  problems <- character(0)
  for (nm in names(pool)) {
    A <- adjacency(pool[[nm]])
    n <- nrow(A)
    for (k in c(0, 2)) {
      for (x in list(sample(0:1, n, replace = TRUE), sample(0:4, n, replace = TRUE))) {
        if (!same_profile(A, x, k)) problems <- c(problems, nm)
      }
    }
  }
  expect_identical(problems, character(0))
})

test_that("profiles agree with brute force on larger random graphs", {
  set.seed(16)
  for (i in 1:20) {
    repeat {
      g <- igraph::sample_gnp(sample(8:12, 1), runif(1, 0.15, 0.45))
      if (igraph::is_connected(g)) break
    }
    A <- adjacency(g)
    for (k in c(0, 3)) expect_true(same_profile(A, sample(0:3, nrow(A), replace = TRUE), k))
  }
})

test_that("very large counts are rescaled, and profiles that would lose rare counts are not computed", {
  # The first member of each pair has value 1: among the 2^400 maximum sets,
  # choose(400, r) have value sum r.
  p <- .rkn_family_profile(adjacency(windmill(400)), c(0, rep(c(1, 0), 400)), 0)
  expect_equal(p$alpha, 400)
  expect_equal(log(p$counts[1, ]) + p$log_scale, lchoose(400, 0:400), tolerance = 1e-10)
  # With 600 pairs the counts span more than 300 natural orders of magnitude.
  expect_null(.rkn_family_profile(adjacency(windmill(600)), c(0, rep(c(1, 0), 600)), 0))
})

test_that("the budget is checked before the profile is computed", {
  g <- igraph::make_lattice(c(4, 5))
  A <- adjacency(g)
  x <- rep(0:2, length.out = nrow(A))
  elim <- .rkn_elimination_order(A)
  budget <- .rkn_profile_budget(elim, x, 2)
  expect_false(is.null(.rkn_family_profile(A, x, 2, elim, max_cells = budget[["cells"]], max_work = budget[["work"]])))
  expect_null(.rkn_family_profile(A, x, 2, elim, max_cells = budget[["cells"]] - 1))
  expect_null(.rkn_family_profile(A, x, 2, elim, max_work = budget[["work"]] - 1))
})

test_that("partial profiles find the smallest and largest value sums among the maximum sets", {
  width_limit_only()
  set.seed(17)
  for (i in 1:20) {
    repeat {
      g <- igraph::sample_gnp(sample(8:12, 1), runif(1, 0.15, 0.45))
      if (igraph::is_connected(g)) break
    }
    A <- adjacency(g)
    x <- sample(0:5, nrow(A), replace = TRUE)
    b <- profile_brute(A, x, 0)
    reachable <- which(b$counts[1, ] > 0) - 1
    p <- .rkn_partial_profile(A, x, n_draws = 20)
    expect_equal(p$alpha, b$alpha)
    expect_equal(range(p$values), range(reachable))
    expect_true(all(p$values %in% reachable))
    for (j in seq_along(p$sets)) {
      s <- p$sets[[j]]
      expect_equal(sum(s), b$alpha)
      expect_false(any(A[s, s]))
      expect_equal(sum(x[s]), p$values[j])
    }
  }
  # Too wide for the engine: weighted dominance keeps the extremes exact. The hub (9)
  # is dominated by 6, which has the same value, whichever extreme is sought.
  fan <- igraph::make_ring(8, circular = FALSE) + igraph::vertex(1)
  fan <- igraph::add_edges(fan, as.vector(rbind(9, 1:8)))
  A <- adjacency(fan)
  x <- c(3, 0, 1, 4, 0, 2, 5, 1, 2)
  reachable <- which(profile_brute(A, x, 0)$counts[1, ] > 0) - 1
  expect_equal(range(.rkn_partial_profile(A, x, n_draws = 5, max_width = 1)$values), range(reachable))
  expect_null(.rkn_partial_profile(adjacency(igraph::make_ring(7)), rep(1, 7), max_width = 1))
})

test_that("draws from a profile are uniform among the sets of the given size and value", {
  set.seed(18)
  for (i in 1:6) {
    repeat {
      g <- igraph::sample_gnp(sample(6:9, 1), runif(1, 0.2, 0.45))
      if (igraph::is_connected(g)) break
    }
    A <- adjacency(g)
    n <- nrow(A)
    x <- sample(0:2, n, replace = TRUE)
    prof <- .rkn_family_profile(A, x, 2, keep = TRUE)
    # The most frequent combination of size and value, and its sets.
    labels <- character(0)
    keys <- character(0)
    for (s in 0:(2^n - 1)) {
      sel <- bitwAnd(s, 2^(seq_len(n) - 1)) > 0
      if (sum(sel) > 1 && any(A[sel, sel])) next
      if (prof$alpha - sum(sel) > 2) next
      keys <- c(keys, paste(prof$alpha - sum(sel), sum(x[sel])))
      labels <- c(labels, paste(which(sel), collapse = ","))
    }
    key <- names(which.max(table(keys)))
    target <- labels[keys == key]
    jr <- as.numeric(strsplit(key, " ")[[1]])
    drawn <- character(150 * length(target))
    for (d in seq_along(drawn)) drawn[d] <- paste(which(.rkn_profile_sample(prof, jr[1], jr[2])), collapse = ",")
    expect_setequal(unique(drawn), target)
    if (length(target) > 1) expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
  }
})

test_that("partial profiles fall back to the reduction without weights", {
  width_limit_only()
  # The centre has the largest value, so with the weights of the largest sum no leaf
  # dominates it; without them the leaves do, and the only maximum set is left.
  star <- igraph::make_star(4, mode = "undirected")
  part <- .rkn_partial_profile(adjacency(star), c(5, 0, 1, 2), n_draws = 3, max_width = 0)
  expect_false(is.null(part))
  expect_identical(part$alpha, 3)
  expect_identical(part$values, 3)
  expect_identical(part$sets[[1]], c(FALSE, TRUE, TRUE, TRUE))
})
