crits <- list(
  target_mean = rkn_criterion_target_mean("age", target = 50, tolerance = 30),
  cutoff = rkn_criterion_cutoff("age", cutoff = 50, prefer = "above"),
  level = rkn_criterion_level("sex", "F"),
  balance = rkn_criterion_category_balance("grp", c(X = 0.5, Y = 0.5))
)

test_that("per-family criteria selection (<= 5 nodes) picks the best maximum set compatible with the a priori decisions", {
  graphs <- atlas_connected(1, 5)
  problems <- character(0)
  for (nm in names(graphs)) {
    g <- graphs[[nm]]
    ref <- brute_force_mis(g)
    for (cov_seed in 1:3) {
      fam_df <- add_covariates(fam_df_from_graph(g), cov_seed)
      for (cn in names(crits)) {
        out <- rkn_select_criteria_upto5(fam_df, crits[[cn]], seed = 1,
                                         return_diagnostics = TRUE, g_sub = g)
        p <- selection_problems(g, out$kept, out$discarded, ref$alpha)
        best <- -Inf
        for (s in ref$sets) {
          if (all(out$kept_apriori %in% s) && !any(out$discarded_apriori %in% s)) {
            best <- max(best, crits[[cn]]$score_fn(fam_df, s))
          }
        }
        got <- crits[[cn]]$score_fn(fam_df, out$kept)
        if (best == -Inf) {
          p <- c(p, "no maximum set is compatible with the a priori decisions")
        } else if (got < best - 1e-12) {
          p <- c(p, sprintf("score %.3f below the best compatible score %.3f", got, best))
        }
        if (length(p) > 0) problems <- c(problems, paste0(nm, " cov ", cov_seed, " ", cn, ": ", p))
      }
    }
  }
  expect_identical(problems, character(0))
})

test_that("per-family criteria selection (6-7 nodes) returns a maximum independent set", {
  graphs <- thin_on_cran(atlas_connected(6, 7))
  problems <- character(0)
  for (nm in names(graphs)) {
    g <- graphs[[nm]]
    alpha <- brute_force_mis(g)$alpha
    fam_df <- add_covariates(fam_df_from_graph(g), 1)
    for (cn in names(crits)) {
      out <- rkn_select_criteria_6plus(fam_df, g, crits[[cn]], seed = 1)
      p <- selection_problems(g, out$kept, out$discarded, alpha)
      if (length(p) > 0) problems <- c(problems, paste0(nm, " ", cn, ": ", p))
    }
  }
  expect_identical(problems, character(0))
})

test_that("on a 6-clique the criterion keeps the preferred node", {
  g <- igraph::make_full_graph(6)
  igraph::V(g)$name <- as.character(1:6)
  fam_df <- fam_df_from_graph(g)
  fam_df$sex <- c("M", "M", "M", "F", "M", "M")
  out <- rkn_select_criteria_6plus(fam_df, g, rkn_criterion_level("sex", "F"), seed = 1)
  expect_identical(out$kept, "4")
})

test_that("on an even path the criterion keeps the better of the two alternating sides", {
  g <- igraph::make_ring(6, circular = FALSE)
  igraph::V(g)$name <- as.character(1:6)
  fam_df <- fam_df_from_graph(g)
  fam_df$sex <- c("M", "F", "M", "M", "M", "M")
  out <- rkn_select_criteria_6plus(fam_df, g, rkn_criterion_level("sex", "F"), seed = 1)
  expect_setequal(out$kept, c("2", "4", "6"))
})

test_that("a secondary criterion breaks ties left by the primary one", {
  g <- igraph::make_full_graph(3)
  igraph::V(g)$name <- c("a", "b", "c")
  fam_df <- fam_df_from_graph(g)
  fam_df$sex <- c("F", "F", "M")
  fam_df$age <- c(30, 68, 70)
  out <- rkn_select_criteria_upto5(fam_df,
                                   criterion_primary = rkn_criterion_level("sex", "F"),
                                   criterion_secondary = rkn_criterion_target_mean("age", 70, 50),
                                   seed = 1, g_sub = g)
  expect_identical(out$kept, "b")
})

test_that("criteria selectors validate their criterion arguments", {
  g <- igraph::make_full_graph(3)
  igraph::V(g)$name <- c("a", "b", "c")
  fam_df <- fam_df_from_graph(g)
  expect_error(rkn_select_criteria_upto5(fam_df), "criterion_primary")
  expect_error(rkn_select_criteria_upto5(fam_df, crits$level, criterion_secondary = "x"), "criterion_secondary")
  expect_error(rkn_select_criteria(make_network(1), criterion_primary = "x"), "criterion_primary")
})

test_that("network-level criteria selection reads covariates from vertex attributes and returns a maximum independent set", {
  skip_if_not_installed("lpSolve")
  for (seed in 1:3) {
    g <- make_network(seed)
    set.seed(seed)
    igraph::V(g)$age <- sample(20:80, igraph::vcount(g), replace = TRUE)
    igraph::V(g)$sex <- sample(c("F", "M"), igraph::vcount(g), replace = TRUE)
    res <- rkn_select_criteria(g, criterion_primary = crits$level, seed = seed)
    expect_identical(selection_problems(g, res$kept, res$discarded, network_alpha(g)), character(0))
  }
})

test_that("network-level criteria selection is reproducible", {
  g <- make_network(21)
  set.seed(21)
  igraph::V(g)$sex <- sample(c("F", "M"), igraph::vcount(g), replace = TRUE)
  expect_identical(rkn_select_criteria(g, criterion_primary = crits$level, seed = 2),
                   rkn_select_criteria(g, criterion_primary = crits$level, seed = 2))
})

test_that("the criteria selection records its criteria in the rkn_selection object", {
  g <- make_network(22)
  set.seed(22)
  igraph::V(g)$sex <- sample(c("F", "M"), igraph::vcount(g), replace = TRUE)
  igraph::V(g)$age <- sample(20:80, igraph::vcount(g), replace = TRUE)
  sel <- rkn_select_criteria(g, crits$level, crits$target_mean, seed = 2, diagnostics = FALSE)
  expect_s3_class(sel, "rkn_selection")
  expect_identical(sel$method, "criteria")
  expect_identical(sel$criteria, list(primary = crits$level, secondary = crits$target_mean))
  expect_identical(sel$seed, 2)
  expect_false(sel$diagnostics)
})
