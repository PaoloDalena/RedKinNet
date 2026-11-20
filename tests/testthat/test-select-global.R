best_global <- function(g, k, p) {
  sex <- igraph::V(g)$sex
  miss <- as.numeric(is.na(sex))
  cat <- as.numeric(!is.na(sex) & sex == "F")
  ids <- igraph::V(g)$name
  best <- best_selections(g, cat, miss, k, p)
  labels <- character(length(best))
  for (i in seq_along(best)) {
    labels[i] <- paste(sort(ids[as.integer(strsplit(best[i], ",")[[1]])]), collapse = ",")
  }
  labels
}

test_that("global selections reach the best proportion of the whole selection, with the maximum size", {
  g <- example_network()
  sel <- rkn_select_global(g, "sex", c(F = 0.5), seed = 1)
  expect_s3_class(sel, "rkn_selection")
  expect_identical(sel$method, "global")
  best <- best_global(g, 0, 0.5)
  expect_true(paste(sort(sel$kept), collapse = ",") %in% best)
  expect_equal(exp(sel$global$log_count), length(best))
  expect_identical(sel$global$size, 8)
  expect_identical(sel$global$n_category, 4)

  set.seed(22)
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1))
  for (i in 1:15) {
    repeat {
      parts <- pool[sample(length(pool), sample(2:4, 1))]
      for (j in seq_along(parts)) igraph::V(parts[[j]])$name <- paste0(j, "_", seq_len(igraph::vcount(parts[[j]])))
      net <- do.call(igraph::disjoint_union, unname(parts))
      if (igraph::vcount(net) <= 12) break
    }
    igraph::V(net)$sex <- sample(c("F", "M", NA), igraph::vcount(net), replace = TRUE, prob = c(0.45, 0.45, 0.1))
    k <- sample(0:2, 1)
    p <- stats::runif(1)
    sel <- rkn_select_global(net, "sex", c(F = p), max_loss = k, seed = i)
    best <- best_global(net, k, p)
    expect_true(paste(sort(sel$kept), collapse = ",") %in% best)
    expect_equal(exp(sel$global$log_count), length(best))
  }
})

test_that("global selections are reproducible and leave the caller's random number stream untouched", {
  g <- example_network()
  expect_identical(rkn_select_global(g, "sex", c(F = 0.5), seed = 3)$kept,
                   rkn_select_global(g, "sex", c(F = 0.5), seed = 3)$kept)
  set.seed(123)
  rkn_select_global(g, "sex", c(F = 0.5), seed = 3)
  after <- stats::runif(3)
  set.seed(123)
  expect_identical(after, stats::runif(3))
})

test_that("the target defaults to the proportion of the network", {
  g <- example_network()
  sel <- rkn_select_global(g, "sex", seed = 1)
  expect_true(sel$global$cohort)
  expect_equal(sel$global$target, 7 / 15)
  expect_identical(sel$global$category, "F")
  expect_output(print(sel), "target 0.4667, the proportion of the network", fixed = TRUE)
})

test_that("the proportion is computed among the selected individuals with a known value", {
  g <- example_network()
  igraph::V(g)$sex[igraph::V(g)$name %in% c("U1", "P2")] <- NA
  sel <- rkn_select_global(g, "sex", c(F = 0.5), seed = 1)
  expect_identical(sel$global$n_missing, 2)
  expect_identical(sel$global$n_known, 6)
  expect_equal(sel$global$proportion, sel$global$n_category / 6)
  expect_output(print(sel), "with a known value")
  expect_output(print(sel), "Selected with a missing value of sex: 2")
  expect_identical(sum(sel$family_summary$n_missing), 2L)
})

test_that("max_loss gives up individuals to bring the proportion closer to the target", {
  g <- example_network()
  full <- rkn_select_global(g, "sex", c(F = 0.7), seed = 1)
  expect_identical(full$global$size, 8)
  expect_equal(full$global$proportion, 5 / 8)
  lossy <- rkn_select_global(g, "sex", c(F = 0.7), max_loss = 2, seed = 1)
  expect_identical(lossy$global$size, 7)
  expect_identical(lossy$global$max_size, 8)
  expect_equal(lossy$global$proportion, 5 / 7)
  expect_equal(exp(lossy$global$log_count), length(best_global(g, 2, 0.7)))
  expect_output(print(lossy), "Selected: 7 of a largest possible 8 (max_loss = 2)", fixed = TRUE)
  # No loss when the target is already reached at the maximum size.
  expect_identical(rkn_select_global(g, "sex", c(F = 0.5), max_loss = 2, seed = 1)$global$size, 8)
  # A max_loss beyond the size of the network changes nothing and keeps every family exact.
  huge <- rkn_select_global(g, "sex", c(F = 0.7), max_loss = 1e6, seed = 1)
  whole <- rkn_select_global(g, "sex", c(F = 0.7), max_loss = 15, seed = 1)
  expect_true(huge$global$exact)
  expect_identical(huge$global$n_inexact, 0L)
  expect_identical(huge$kept, whole$kept)
  expect_equal(huge$global$log_count, whole$global$log_count)
  expect_identical(huge$global$max_loss, 1e6)
})

test_that("a message warns before a selection beyond the computational budget", {
  real <- .rkn_global_select
  local_mocked_bindings(.rkn_global_select = function(...) real(..., max_work = 1))
  g <- example_network()
  expect_message(sel <- rkn_select_global(g, "sex", c(F = 0.5), seed = 1), "may take several minutes")
  expect_identical(sel$global$n_category, 4)
})

test_that("the inputs of a global selection are validated", {
  g <- example_network()
  expect_error(rkn_select_global(g), "covariate")
  expect_error(rkn_select_global(g, "sex", c(F = 0.5), max_loss = -1), "max_loss")
  expect_error(rkn_select_global(g, "sex", c(F = 0.5), max_loss = 1.5), "max_loss")
  expect_error(rkn_select_global(g, "age", c(F = 0.5)), "no vertex attribute")
  expect_error(rkn_select_global(g, "sex", c(X = 0.5)), "not found")
  expect_error(rkn_select_global(g, "sex", 0.5), "name the category")
  expect_error(rkn_select_global(g, "sex", c(F = 1.2)), "between 0 and 1")
  expect_error(rkn_select_global(g, "sex", c(F = 0.3, M = 0.3)), "sum to 1")
  igraph::V(g)$age <- seq_len(igraph::vcount(g))
  expect_error(rkn_select_global(g, "age", c(F = 0.5)), "single unnamed number")
  expect_error(rkn_select_global(g, "age", c(40, 50)), "single unnamed number")
  expect_error(rkn_select_global(g, "age", 10, resolution = 0), "single positive number")
  expect_error(rkn_select_global(g, "sex", c(F = 0.5), resolution = 1), "numeric covariates only")
  igraph::V(g)$age[1] <- Inf
  expect_error(rkn_select_global(g, "age", 10), "infinite")
  igraph::V(g)$sex[1:3] <- "other"
  expect_error(rkn_select_global(g, "sex"), "more than two categories")
  expect_error(rkn_select_global(g, "sex", c(F = 0.5, M = 0.5)), "one category")
  expect_s3_class(rkn_select_global(g, "sex", c(F = 0.5), seed = 1), "rkn_selection")
})

test_that("factor and logical covariates are categories", {
  g <- example_network()
  sel <- rkn_select_global(g, "sex", c(F = 0.5), seed = 2)
  igraph::V(g)$sex <- factor(igraph::V(g)$sex)
  expect_identical(rkn_select_global(g, "sex", c(F = 0.5), seed = 2)$kept, sel$kept)
  igraph::V(g)$woman <- as.character(igraph::V(g)$sex) == "F"
  expect_identical(rkn_select_global(g, "woman", c(`TRUE` = 0.5), seed = 2)$kept, sel$kept)
  # The default target is the proportion of TRUE, or of the first level of a factor.
  expect_identical(rkn_select_global(g, "woman", seed = 1)$global$category, "TRUE")
  igraph::V(g)$sex <- factor(igraph::V(g)$sex, levels = c("M", "F"))
  expect_identical(rkn_select_global(g, "sex", seed = 1)$global$category, "M")
  expect_equal(rkn_select_global(g, "sex", seed = 1)$global$target, 8 / 15)
})

test_that("global selections store family sizes and identifiers as integers, as the others do", {
  sel <- rkn_select_global(example_network(), "sex", c(F = 0.5), seed = 1)
  expect_type(sel$family_sizes, "integer")
  expect_type(sel$family_summary$family_size, "integer")
  expect_type(sel$node_diagnostics$family_id, "integer")
  expect_type(sel$node_diagnostics$family_size, "integer")
})

test_that("print, summary and plot describe global selections", {
  g <- example_network()
  sel <- rkn_select_global(g, "sex", c(F = 0.5), seed = 1)
  expect_output(print(sel), "Global selection (proportion of sex = F), seed 1", fixed = TRUE)
  expect_output(print(sel), "Proportion of sex = F among the selected: 0.5000 (4 of 8; target 0.5)", fixed = TRUE)
  expect_output(print(sel), "Selections reaching this proportion: 3, one drawn uniformly", fixed = TRUE)
  expect_silent(capture.output(print(summary(sel))))
  expect_named(sel$family_summary, c("family_id", "family_size", "selector", "exact", "n_kept", "n_discarded",
                                     "n_category", "n_missing"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_message(plot(sel, type = "diagnostics"), "not available for global selections")

  inexact <- sel
  inexact$global$exact <- FALSE
  inexact$global$n_inexact <- 2
  expect_output(print(inexact), "Families beyond the computational budget: 2")

  plain <- rkn_select_global(g, "sex", c(F = 0.5), seed = 1, diagnostics = FALSE)
  expect_null(plain$family_summary)
  expect_identical(plain$kept, sel$kept)
  expect_message(expect_output(print(plain), "Proportion of sex = F"), "diagnostics = FALSE")
})

labels_of <- function(best, ids) {
  out <- character(length(best))
  for (i in seq_along(best)) out[i] <- paste(sort(ids[as.integer(strsplit(best[i], ",")[[1]])]), collapse = ",")
  out
}

random_small_network <- function(pool) {
  repeat {
    parts <- pool[sample(length(pool), sample(2:4, 1))]
    for (j in seq_along(parts)) igraph::V(parts[[j]])$name <- paste0(j, "_", seq_len(igraph::vcount(parts[[j]])))
    net <- do.call(igraph::disjoint_union, unname(parts))
    if (igraph::vcount(net) <= 12) return(net)
  }
}

test_that("global selections reach the best mean of a numeric covariate", {
  set.seed(23)
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1))
  for (i in 1:15) {
    net <- random_small_network(pool)
    age <- sample(c(20:70, NA), igraph::vcount(net), replace = TRUE, prob = c(rep(0.9 / 51, 51), 0.1))
    igraph::V(net)$age <- age
    k <- sample(0:2, 1)
    p <- stats::runif(1, 20, 70)
    sel <- rkn_select_global(net, "age", p, max_loss = k, seed = i)
    expect_false(sel$global$rounded)
    best <- labels_of(best_selections(net, ifelse(is.na(age), 0, age), as.numeric(is.na(age)), k, p),
                      igraph::V(net)$name)
    expect_true(paste(sort(sel$kept), collapse = ",") %in% best)
    expect_equal(exp(sel$global$log_count), length(best))
    expect_equal(sel$global$mean, mean(age[igraph::V(net)$name %in% sel$kept], na.rm = TRUE))
  }
})

test_that("rounded values give the best selection for them, within one resolution of the best", {
  set.seed(24)
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1))
  for (i in 1:10) {
    net <- random_small_network(pool)
    bmi <- round(stats::runif(igraph::vcount(net), 17, 40), 1)
    igraph::V(net)$bmi <- bmi
    p <- stats::runif(1, 20, 35)
    sel <- rkn_select_global(net, "bmi", p, max_loss = 1, resolution = 2, seed = i)
    expect_true(sel$global$rounded)
    steps <- round((bmi - min(bmi)) / 2)
    best <- labels_of(best_selections(net, steps, rep(0, length(bmi)), 1, (p - min(bmi)) / 2), igraph::V(net)$name)
    expect_true(paste(sort(sel$kept), collapse = ",") %in% best)
    exact <- labels_of(best_selections(net, bmi, rep(0, length(bmi)), 1, p), igraph::V(net)$name)[1]
    best_mean <- mean(bmi[igraph::V(net)$name %in% strsplit(exact, ",")[[1]]])
    expect_lte(abs(sel$global$mean - p), abs(best_mean - p) + 2)
  }
})

test_that("the default resolution keeps integers and divides other ranges into at most 100 steps", {
  expect_identical(.rkn_default_resolution(20:85), 1)
  expect_identical(.rkn_default_resolution(c(0, 1)), 1)
  expect_identical(.rkn_default_resolution(c(5, 5)), 1)
  expect_equal(.rkn_default_resolution(c(16.2, 45.1)), 0.5)
  expect_equal(.rkn_default_resolution(c(0, 3.3)), 0.05)
  expect_equal(.rkn_default_resolution(c(0, 150)), 2)
  expect_equal(.rkn_default_resolution(c(0, 200000)), 2000)
})

test_that("print and the family summary describe selections by mean", {
  g <- example_network()
  igraph::V(g)$age <- c(52, 21, 55, 34, 31, 60, 33, 29, 63, 35, 61, 38, 30, 58, 44)[
    match(igraph::V(g)$name, c("P1", "C1", "P2", "S1", "S2", "T1", "T2", "T3", "D1", "K1", "D2", "K2", "K3", "D3", "U1"))]
  sel <- rkn_select_global(g, "age", 45, seed = 1)
  expect_identical(sel$global$type, "mean")
  expect_equal(sel$global$mean, 393 / 8)
  expect_output(print(sel), "Global selection (mean of age), seed 1", fixed = TRUE)
  expect_output(print(sel), "Mean of age among the selected: 49.12 (8 individuals; target 45)", fixed = TRUE)
  expect_named(sel$family_summary, c("family_id", "family_size", "selector", "exact", "n_kept", "n_discarded",
                                     "mean", "n_missing"))
  expect_equal(sel$family_summary$mean[1], 53.5)
  expect_output(print(rkn_select_global(g, "age", seed = 1)), "the mean of the network")

  igraph::V(g)$bmi <- c(22.3, 31.7, 27.4, 19.8, 24.1, 29.9, 21.5, 33.2, 26.6, 23.4, 28.8, 25.0, 20.7, 30.4, 24.9)
  rounded <- rkn_select_global(g, "bmi", 25, seed = 1)
  expect_equal(rounded$global$resolution, 0.2)
  expect_output(print(rounded), "Values of bmi rounded to multiples of 0.2 for the search")
  expect_output(print(rounded), "Selections reaching this mean of the rounded values: 1")
})

test_that("an error explains when no selection includes an individual with a known value", {
  # A star whose centre alone has a known value: the largest selections are the leaves.
  g <- igraph::make_star(4, mode = "undirected")
  igraph::V(g)$name <- c("hub", "a", "b", "c")
  igraph::V(g)$sex <- c("F", NA, NA, NA)
  expect_error(rkn_select_global(g, "sex", c(F = 0.5)),
               "No selection of the largest size includes an individual with a known value of \"sex\"")
  expect_error(rkn_select_global(g, "sex", c(F = 0.5), max_loss = 1), "within `max_loss`")
  expect_identical(rkn_select_global(g, "sex", c(F = 0.5), max_loss = 2, seed = 1)$kept, "hub")
  expect_error(rkn_tradeoff(g, "sex", c(F = 0.5), max_loss = 1),
               "giving up at most 1 individual includes an individual with a known value of \"sex\"")
})

test_that("the default resolution is coarsened until the tables fit in the memory budget", {
  # 30 related pairs with ages from 1 to 99: the tables need 15,934 cells at resolution 1,
  # 7,973 at 2 and 3,194 at 5, so a budget of 5,000 cells leads to a resolution of 5.
  g <- igraph::make_empty_graph(60, directed = FALSE)
  g <- igraph::add_edges(g, as.vector(rbind(seq(1, 59, 2), seq(2, 60, 2))))
  igraph::V(g)$name <- paste0("s", 1:60)
  igraph::V(g)$age <- as.numeric(c(rbind(1:30, 100 - 1:30)))
  real <- .rkn_global_select
  local_mocked_bindings(.rkn_global_select = function(fams, k, ...) real(fams, k, ..., max_cells = 5000))
  sel <- rkn_select_global(g, "age", 50, seed = 1)
  expect_identical(sel$global$resolution, 5)
  expect_true(sel$global$coarsened)
  expect_length(sel$kept, 30)
  expect_output(print(sel), "coarser than the default to fit the memory budget")
  expect_error(rkn_select_global(g, "age", 50, resolution = 1, seed = 1), class = "rkn_table_limit")
  expect_false(rkn_select_global(g, "age", 50, resolution = 5, seed = 1)$global$coarsened)
  # The trade-off keeps fewer tables at once (about 4,100 cells at resolution 1).
  expect_false(attr(rkn_tradeoff(g, "age", 50, max_loss = 0), "tradeoff")$coarsened)
  local_mocked_bindings(.rkn_global_select = function(fams, k, ...) real(fams, k, ..., max_cells = 2000))
  tr <- rkn_tradeoff(g, "age", 50, max_loss = 0)
  expect_identical(attr(tr, "tradeoff")$resolution, 5)
  expect_true(attr(tr, "tradeoff")$coarsened)
  expect_output(print(tr), "coarser than the default to fit the memory budget")
})

test_that("coarser resolutions follow the sequence 1, 2, 5 and exceed the previous one", {
  expect_identical(.rkn_coarser_resolution(1, 3.19), 5)
  expect_identical(.rkn_coarser_resolution(1, 6.4), 10)
  expect_identical(.rkn_coarser_resolution(0.5, 1.1), 1)
  expect_identical(.rkn_coarser_resolution(2, 1.0001), 5)
})
