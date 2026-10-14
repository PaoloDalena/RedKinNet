# Until the selectors return rkn_selection objects, test objects are assembled from
# their current outputs.
make_selection <- function(diagnostics = TRUE, method = "random", criteria = NULL) {
  g <- make_network(51)
  res <- rkn_select_random(g, seed = 1, diagnostics = TRUE, verbose = FALSE)
  .new_rkn_selection(
    kept = res$kept,
    discarded = res$discarded,
    family_sizes = as.integer(igraph::components(g)$csize),
    method = method,
    criteria = criteria,
    seed = 1,
    node_diagnostics = if (diagnostics) res$node_diagnostics else NULL,
    family_summary = if (diagnostics) res$family_summary else NULL,
    graph = g,
    call = quote(rkn_select_random(g = g, seed = 1))
  )
}

test_that("the constructor returns an object of class rkn_selection", {
  sel <- make_selection()
  expect_s3_class(sel, "rkn_selection")
  expect_true(sel$diagnostics)
  expect_false(make_selection(diagnostics = FALSE)$diagnostics)
})

test_that("print shows the method, the call and the selection counts, and returns x invisibly", {
  sel <- make_selection()
  expect_output(print(sel), "<rkn_selection> Random selection, seed 1")
  expect_output(print(sel), "Call: rkn_select_random(g = g, seed = 1)", fixed = TRUE)
  expect_output(print(sel), sprintf("selected: %d", length(sel$kept)))
  expect_output(print(sel), "enumerating maximum independent sets")
  capture.output(out <- withVisible(print(sel)))
  expect_false(out$visible)
  expect_identical(out$value, sel)
})

test_that("print names the criteria of a criteria-based selection", {
  crit <- list(primary = rkn_criterion_level("sex", "F"),
               secondary = rkn_criterion_target_mean("age", 50, 10))
  sel <- make_selection(method = "criteria", criteria = crit)
  expect_output(print(sel), "Criteria-based selection (level, then target_mean)", fixed = TRUE)
})

test_that("print signals when diagnostics were not computed", {
  expect_no_message(expect_output(print(make_selection())))
  expect_message(expect_output(print(make_selection(diagnostics = FALSE))), "diagnostics = TRUE")
})
