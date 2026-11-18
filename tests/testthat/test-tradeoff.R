# For each number of individuals given up, the range of proportions of women, the
# closest to p and the number of selections reaching it, by exhaustive enumeration.
brute_tradeoff <- function(g, p, x = NULL) {
  sex <- igraph::V(g)$sex
  miss <- is.na(sex)
  cat <- !miss & sex == "F"
  if (!is.null(x)) {
    miss <- is.na(x)
    cat <- ifelse(miss, 0, x)
  }
  A <- .rkn_adjacency(g)
  n <- nrow(A)
  size <- C <- M <- numeric(0)
  for (s in 0:(2^n - 1)) {
    sel <- bitwAnd(s, 2^(seq_len(n) - 1)) > 0
    if (sum(sel) > 1 && any(A[sel, sel])) next
    size <- c(size, sum(sel))
    C <- c(C, sum(cat[sel]))
    M <- c(M, sum(miss[sel]))
  }
  alpha <- max(size)
  out <- data.frame(loss = 0:alpha, proportion = NA_real_, distance = NA_real_, min_proportion = NA_real_,
                    max_proportion = NA_real_, count = NA_real_)
  for (j in 0:alpha) {
    at <- size == alpha - j & size - M > 0
    if (!any(at)) next
    prop <- C[at] / (size[at] - M[at])
    d <- abs(prop - p)
    best <- d <= min(d) + 1e-12
    out[j + 1, -1] <- c(min(prop[best]), min(d), min(prop), max(prop), sum(best))
  }
  out
}

singletons <- function(n_women, n_men) {
  g <- igraph::make_empty_graph(n_women + n_men, directed = FALSE)
  igraph::V(g)$name <- paste0("s", seq_len(n_women + n_men))
  igraph::V(g)$sex <- rep(c("F", "M"), c(n_women, n_men))
  g
}

test_that("the trade-off table matches exhaustive enumeration", {
  check <- function(g, p) {
    ref <- brute_tradeoff(g, p)
    tr <- rkn_tradeoff(g, "sex", c(F = p), max_loss = nrow(ref) - 1)
    expect_s3_class(tr, "rkn_tradeoff")
    expect_identical(tr$loss, 0:(nrow(ref) - 1))
    expect_equal(tr$size, as.integer(rev(tr$loss)))
    for (col in c("proportion", "distance", "min_proportion", "max_proportion")) expect_equal(tr[[col]], ref[[col]])
    expect_equal(exp(tr$log_count), ref$count)
  }
  check(example_network(), 0.7)

  set.seed(31)
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1))
  for (i in 1:12) {
    repeat {
      parts <- pool[sample(length(pool), sample(2:4, 1))]
      for (j in seq_along(parts)) igraph::V(parts[[j]])$name <- paste0(j, "_", seq_len(igraph::vcount(parts[[j]])))
      net <- do.call(igraph::disjoint_union, unname(parts))
      if (igraph::vcount(net) <= 12) break
    }
    igraph::V(net)$sex <- sample(c("F", "M", NA), igraph::vcount(net), replace = TRUE, prob = c(0.45, 0.45, 0.1))
    check(net, stats::runif(1))
  }
})

test_that("by default the table stops where the target comes within reach", {
  g <- example_network()
  tr <- rkn_tradeoff(g, "sex", c(F = 0.7))
  expect_identical(tr$loss, 0:1)
  expect_identical(attr(tr, "tradeoff")$reached, 1L)
  expect_identical(attr(tr, "tradeoff")$max_size, 8)

  cohort <- rkn_tradeoff(g, "sex")
  expect_identical(nrow(cohort), 1L)
  expect_identical(attr(cohort, "tradeoff")$reached, 0L)

  # 5% of a largest selection of 8, rounded up, is 1.
  far <- rkn_tradeoff(g, "sex", c(F = 1))
  expect_identical(far$loss, 0:1)
  expect_true(is.na(attr(far, "tradeoff")$reached))

  # Beyond the first 10 individuals given up, the table is extended: 108 women and 120
  # men reach half women by giving up 12 men, 5% of 228 rounded up.
  tr <- rkn_tradeoff(singletons(108, 120), "sex", c(F = 0.5))
  expect_identical(attr(tr, "tradeoff")$reached, 12L)
  expect_identical(nrow(tr), 13L)
  given <- rkn_tradeoff(singletons(108, 120), "sex", c(F = 0.5), max_loss = 12)
  for (col in names(tr)) expect_equal(tr[[col]], given[[col]])
  # The same with families, whose options are reused while the table is extended.
  net <- make_network(7, n_atlas = 30, tree_sizes = c(9, 13, 20))
  igraph::V(net)$sex <- rep(c("F", "M", "M"), length.out = igraph::vcount(net))
  steps <- rkn_tradeoff(net, "sex", c(F = 0.5))
  reached <- attr(steps, "tradeoff")$reached
  expect_gt(reached, 1)
  given <- rkn_tradeoff(net, "sex", c(F = 0.5), max_loss = reached)
  for (col in names(steps)) expect_equal(steps[[col]], given[[col]])
  # 105 women and 120 men would need 15, beyond 5% of 225 rounded up.
  short <- rkn_tradeoff(singletons(105, 120), "sex", c(F = 0.5))
  expect_identical(short$loss, 0:12)
  expect_true(is.na(attr(short, "tradeoff")$reached))
  expect_output(print(short), "Target not within reach by giving up up to 12 (5% of the largest selection, rounded up)",
                fixed = TRUE)
})

test_that("a given max_loss covers every number of individuals given up, up to the largest selection", {
  g <- example_network()
  expect_identical(rkn_tradeoff(g, "sex", c(F = 0.7), max_loss = 4)$loss, 0:4)
  expect_identical(rkn_tradeoff(g, "sex", c(F = 0.5), max_loss = 0)$loss, 0L)
  whole <- rkn_tradeoff(g, "sex", c(F = 0.7), max_loss = 100)
  expect_identical(whole$loss, 0:8)
  expect_true(is.na(whole$proportion[9]))
})

test_that("the table predicts the global selection for each max_loss", {
  set.seed(4)
  for (i in 1:6) {
    g <- make_network(i, n_atlas = 8, tree_sizes = c(9, 13))
    igraph::V(g)$sex <- sample(c("F", "M", NA), igraph::vcount(g), replace = TRUE, prob = c(0.4, 0.55, 0.05))
    p <- stats::runif(1)
    tr <- rkn_tradeoff(g, "sex", c(F = p), max_loss = 4)
    for (L in 0:4) {
      sel <- rkn_select_global(g, "sex", c(F = p), max_loss = L, seed = 1)
      rows <- tr[seq_len(L + 1), ]
      d <- min(rows$distance)
      expect_equal(abs(sel$global$proportion - p), d)
      expect_identical(as.integer(sel$global$size), max(rows$size[rows$distance <= d + 1e-12]))
    }
  }
})

test_that("the table stops before exceeding the memory limit when max_loss is NULL", {
  real <- .rkn_global_select
  local_mocked_bindings(.rkn_global_select = function(fams, k, ...) {
    if (k > 8) stop(structure(class = c("rkn_table_limit", "error", "condition"),
                               list(message = "Table limit.", call = NULL)))
    real(fams, k, ...)
  })
  g <- singletons(108, 120)
  # The table is extended to 1, 2, 4, 8 and then 12 individuals given up.
  tr <- rkn_tradeoff(g, "sex", c(F = 0.5))
  expect_identical(tr$loss, 0:8)
  expect_identical(attr(tr, "tradeoff")$limited, "memory")
  expect_output(print(tr), "a longer table would exceed the memory limit")
  expect_error(rkn_tradeoff(g, "sex", c(F = 0.5), max_loss = 12), "Table limit", class = "rkn_table_limit")
})

test_that("the table stops before exceeding the computational budget when max_loss is NULL", {
  real <- .rkn_global_select
  local_mocked_bindings(.rkn_global_select = function(fams, k, ...) {
    real(fams, k, ..., max_work = if (k >= 4) 1 else Inf)
  })
  g <- singletons(108, 120)
  tr <- rkn_tradeoff(g, "sex", c(F = 0.5))
  expect_identical(tr$loss, 0:2)
  expect_identical(attr(tr, "tradeoff")$limited, "work")
  expect_output(print(tr), "a longer table would exceed the computational budget: set max_loss to extend it")
  expect_message(given <- rkn_tradeoff(g, "sex", c(F = 0.5), max_loss = 12), "may take several minutes")
  expect_identical(given$loss, 0:12)
})

test_that("the inputs of the trade-off are validated", {
  g <- example_network()
  expect_error(rkn_tradeoff(g), "covariate")
  expect_error(rkn_tradeoff(g, "age"), "no vertex attribute")
  for (bad in list(-1, 1.5, "a", c(1, 2), NA)) {
    expect_error(rkn_tradeoff(g, "sex", c(F = 0.5), max_loss = bad), "max_loss` must be NULL")
  }
  expect_error(rkn_tradeoff(g, "sex", c(X = 0.5)), "not found")
  expect_error(rkn_tradeoff(g, "sex", c(F = 0.5), resolution = 1), "numeric covariates only")
  igraph::V(g)$age <- seq_len(igraph::vcount(g))
  expect_error(rkn_tradeoff(g, "age", c(F = 45)), "single unnamed number")
})

test_that("print and plot describe the trade-off", {
  g <- example_network()
  tr <- rkn_tradeoff(g, "sex", c(F = 0.7), max_loss = 4)
  expect_output(print(tr), "Proportion of sex = F by individuals given up (target 0.7)", fixed = TRUE)
  expect_output(print(tr), "Largest selection: 8 individuals, proportion of F from 0.3750 to 0.6250", fixed = TRUE)
  expect_output(print(tr), "Target within reach by giving up 1 (12.5% of the largest selection), a selection of 7",
                fixed = TRUE)
  expect_output(print(tr), "selections")
  expect_output(print(rkn_tradeoff(g, "sex")), "no individual needs to be given up")
  expect_output(print(rkn_tradeoff(g, "sex", c(F = 1), max_loss = 1)),
                "Target not within reach by giving up up to 1 (12.5% of the largest selection)", fixed = TRUE)

  long <- rkn_tradeoff(singletons(108, 120), "sex", c(F = 0.5), max_loss = 30)
  expect_output(print(long), "(10 of 31 rows; the object is a data frame with all of them)", fixed = TRUE)

  inexact <- tr
  attr(inexact, "tradeoff")$n_inexact <- 2
  inexact$log_count <- NA_real_
  expect_output(print(inexact), "Families beyond the computational budget: 2")

  igraph::V(g)$sex[1:2] <- NA
  expect_output(print(rkn_tradeoff(g, "sex", c(F = 0.7))), "with a known value of sex")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(tr))
  expect_silent(plot(rkn_tradeoff(example_network(), "sex")))
  expect_silent(plot(long, main = "Trade-off"))
})

test_that("the trade-off of a numeric covariate matches exhaustive enumeration", {
  set.seed(32)
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1))
  for (i in 1:8) {
    repeat {
      parts <- pool[sample(length(pool), sample(2:4, 1))]
      for (j in seq_along(parts)) igraph::V(parts[[j]])$name <- paste0(j, "_", seq_len(igraph::vcount(parts[[j]])))
      net <- do.call(igraph::disjoint_union, unname(parts))
      if (igraph::vcount(net) <= 12) break
    }
    age <- sample(c(20:70, NA), igraph::vcount(net), replace = TRUE, prob = c(rep(0.9 / 51, 51), 0.1))
    igraph::V(net)$age <- age
    p <- stats::runif(1, 20, 70)
    ref <- brute_tradeoff(net, p, age)
    tr <- rkn_tradeoff(net, "age", p, max_loss = nrow(ref) - 1)
    expect_named(tr, c("loss", "size", "mean", "distance", "min_mean", "max_mean", "log_count"))
    expect_equal(tr$mean, ref$proportion)
    expect_equal(tr$distance, ref$distance)
    expect_equal(tr$min_mean, ref$min_proportion)
    expect_equal(tr$max_mean, ref$max_proportion)
    expect_equal(exp(tr$log_count), ref$count)
  }
})

test_that("print and plot describe the trade-off of a numeric covariate", {
  g <- example_network()
  igraph::V(g)$age <- c(52, 21, 55, 34, 31, 60, 33, 29, 63, 35, 61, 38, 30, 58, 44)[
    match(igraph::V(g)$name, c("P1", "C1", "P2", "S1", "S2", "T1", "T2", "T3", "D1", "K1", "D2", "K2", "K3", "D3", "U1"))]
  tr <- rkn_tradeoff(g, "age", 45, max_loss = 4)
  expect_output(print(tr), "Mean of age by individuals given up (target 45)", fixed = TRUE)
  expect_output(print(tr), "Largest selection: 8 individuals, mean of age from 49.12 to 53.38", fixed = TRUE)
  expect_output(print(tr), "min_mean")
  expect_identical(attr(rkn_tradeoff(g, "age", 45), "tradeoff")$reached, 1L)
  igraph::V(g)$bmi <- c(22.3, 31.7, 27.4, 19.8, 24.1, 29.9, 21.5, 33.2, 26.6, 23.4, 28.8, 25.0, 20.7, 30.4, 24.9)
  expect_output(print(rkn_tradeoff(g, "bmi", 25, max_loss = 2)),
                "Values of bmi rounded to multiples of 0.2: each mean is within 0.1 of that of the recorded values",
                fixed = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(tr))
})
