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
  expect_error(rkn_select_global(g, "age", seed = 1), "categorical covariates")
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
