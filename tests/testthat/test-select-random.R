pool <- thin_on_cran(atlas_connected(1, 7))

test_that("per-family random selection returns a maximum independent set on every connected graph up to 7 nodes", {
  problems <- character(0)
  for (nm in names(pool)) {
    g <- pool[[nm]]
    fam_df <- fam_df_from_graph(g)
    alpha <- brute_force_mis(g)$alpha
    for (seed in 1:3) {
      out <- if (nrow(fam_df) <= 5) {
        rkn_select_random_upto5(fam_df, seed = seed, g_sub = g)
      } else {
        rkn_select_random_6plus(fam_df, g, seed = seed)
      }
      p <- selection_problems(g, out$kept, out$discarded, alpha)
      if (length(p) > 0) problems <- c(problems, paste0(nm, " seed ", seed, ": ", p))
    }
  }
  expect_identical(problems, character(0))
})

test_that("per-family random selection is reproducible and unaffected by diagnostics", {
  problems <- character(0)
  for (nm in names(pool)) {
    g <- pool[[nm]]
    fam_df <- fam_df_from_graph(g)
    run <- function(diag) {
      if (nrow(fam_df) <= 5) {
        rkn_select_random_upto5(fam_df, seed = 7, g_sub = g, return_diagnostics = diag)
      } else {
        rkn_select_random_6plus(fam_df, g, seed = 7, return_diagnostics = diag)
      }
    }
    a <- run(FALSE)
    b <- run(FALSE)
    d <- run(TRUE)
    if (!identical(a$kept, b$kept)) problems <- c(problems, paste0(nm, ": not reproducible"))
    if (!setequal(a$kept, d$kept)) problems <- c(problems, paste0(nm, ": diagnostics change the selection"))
  }
  expect_identical(problems, character(0))
})

test_that("network-level random selection returns a maximum independent set", {
  skip_if_not_installed("lpSolve")
  for (seed in 1:5) {
    g <- make_network(seed)
    res <- rkn_select_random(g, seed = seed, diagnostics = TRUE, verbose = FALSE)
    expect_identical(selection_problems(g, res$kept, res$discarded, network_alpha(g)), character(0))
  }
})

test_that("network-level random selection is reproducible", {
  g <- make_network(11)
  expect_identical(rkn_select_random(g, seed = 3, verbose = FALSE),
                   rkn_select_random(g, seed = 3, verbose = FALSE))
})

test_that("the families argument accepts NULL, a membership vector or a list", {
  g <- make_network(12)
  mem <- igraph::components(g)$membership
  fam_list <- split(igraph::V(g)$name, mem)
  k_null <- rkn_select_random(g, seed = 1, verbose = FALSE)
  expect_identical(sort(rkn_select_random(g, families = mem, seed = 1, verbose = FALSE)), sort(k_null))
  expect_identical(sort(rkn_select_random(g, families = fam_list, seed = 1, verbose = FALSE)), sort(k_null))
})

test_that("diagnostics have the documented structure and agree with the selection", {
  g <- make_network(13)
  plain <- rkn_select_random(g, seed = 1, verbose = FALSE)
  expect_type(plain, "character")

  res <- rkn_select_random(g, seed = 1, diagnostics = TRUE, return_results_by_family = TRUE, verbose = FALSE)
  expect_named(res, c("kept", "discarded", "node_diagnostics", "family_summary", "results_by_family"))
  expect_setequal(res$kept, plain)
  expect_equal(nrow(res$node_diagnostics), igraph::vcount(g))
  expect_equal(nrow(res$family_summary), igraph::components(g)$no)
  expect_length(res$results_by_family, igraph::components(g)$no)
  nd <- res$node_diagnostics
  expect_setequal(nd$id[nd$is_selected], res$kept)
  expect_equal(sum(res$family_summary$n_kept), length(res$kept))
})

test_that("verbose output, plot fallback and input validation behave as documented", {
  g <- make_network(14)
  expect_output(rkn_select_random(g, seed = 1), "rkn_select_random\\(\\) summary")

  grDevices::pdf(NULL)
  expect_warning(rkn_select_random(g, seed = 1, verbose = FALSE, plot = "diagnostics"),
                 "requires diagnostics=TRUE")
  grDevices::dev.off()

  expect_error(rkn_select_random(list()), "igraph object")
  g_unnamed <- igraph::make_ring(5)
  expect_error(rkn_select_random(g_unnamed), "must be defined")
  expect_error(rkn_select_random(g, families = 1:3), "Unsupported 'families' format")
})
