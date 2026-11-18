test_that("the constructor returns an object of class rkn_selection", {
  sel <- make_selection()
  expect_s3_class(sel, "rkn_selection")
  expect_true(sel$diagnostics)
  expect_false(make_selection(diagnostics = FALSE)$diagnostics)
})

test_that("print shows the method, the call and the selection counts, and returns x invisibly", {
  sel <- make_selection()
  expect_output(print(sel), "<rkn_selection> Random selection, seed 1")
  expect_output(print(sel), "Call: rkn_select_random(g = g, seed = 1, diagnostics = diagnostics)", fixed = TRUE)
  expect_output(print(sel), sprintf("selected: %d", length(sel$kept)))
  expect_output(print(sel), "a search over maximum sets")
  capture.output(out <- withVisible(print(sel)))
  expect_false(out$visible)
  expect_identical(out$value, sel)
})

test_that("print names the criteria of a criteria-based selection", {
  g <- make_network(51)
  set.seed(51)
  igraph::V(g)$sex <- sample(c("F", "M"), igraph::vcount(g), replace = TRUE)
  igraph::V(g)$age <- sample(20:80, igraph::vcount(g), replace = TRUE)
  sel <- rkn_select_criteria(g, rkn_criterion_level("sex", "F"),
                             rkn_criterion_cutoff("age", 50, "below"), seed = 1)
  expect_output(print(sel), "Criteria-based selection (level, then cutoff), seed 1", fixed = TRUE)
})

test_that("print signals when diagnostics were not computed", {
  expect_no_message(expect_output(print(make_selection())))
  expect_message(expect_output(print(make_selection(diagnostics = FALSE))), "diagnostics = TRUE")
})

test_that("summary returns a summary.rkn_selection with consistent tables", {
  sel <- make_selection()
  s <- summary(sel)
  n_ind <- length(sel$kept) + length(sel$discarded)
  expect_s3_class(s, "summary.rkn_selection")
  expect_equal(s$overall$n_individuals, n_ind)
  expect_equal(s$overall$n_selected, length(sel$kept))
  expect_equal(sum(s$by_structure$n_families), length(sel$family_sizes))
  expect_equal(sum(s$by_structure$n_individuals), n_ind)
  expect_equal(sum(s$by_structure$n_selected), length(sel$kept))
  expect_equal(sum(s$decisions), n_ind)
  expect_equal(sum(s$decisions[, "selected"]), length(sel$kept))
})

test_that("summary assigns each family to its structure", {
  parts <- list(igraph::make_empty_graph(1, directed = FALSE), igraph::make_ring(2),
                igraph::make_ring(3, circular = FALSE), igraph::make_ring(3), igraph::make_star(4, mode = "undirected"),
                igraph::make_ring(5), igraph::make_ring(7), igraph::make_ring(9, circular = FALSE))
  g <- igraph::make_empty_graph(0, directed = FALSE)
  for (p in parts) g <- igraph::disjoint_union(g, p)
  igraph::V(g)$name <- paste0("id", seq_len(igraph::vcount(g)))
  s <- summary(rkn_select_random(g, seed = 1, diagnostics = FALSE))
  expect_equal(s$by_structure$family_size, c("1", "2", "3", "3", "4", "5", "7-9"))
  expect_equal(s$by_structure$structure, .rkn_structure_classes)
  expect_equal(s$by_structure$n_families, c(1, 1, 1, 1, 1, 1, 2))
  expect_equal(s$by_structure$n_selected, c(1, 1, 2, 1, 3, 2, 8))
})

test_that("summary without diagnostics keeps the structure table and drops the decision table", {
  s <- summary(make_selection(diagnostics = FALSE))
  expect_equal(sum(s$by_structure$n_individuals), s$overall$n_individuals)
  expect_null(s$decisions)
})

test_that("the summary print shows its sections and signals missing diagnostics", {
  s <- summary(make_selection())
  expect_output(print(s), "<rkn_selection summary> Random selection, seed 1")
  expect_output(print(s), "Selection by family structure")
  expect_output(print(s), "Individuals by type of decision and outcome")
  expect_output(print(s), "decided among alternative sets")
  capture.output(out <- withVisible(print(s)))
  expect_false(out$visible)
  expect_message(expect_output(print(summary(make_selection(diagnostics = FALSE)))), "diagnostics = TRUE")
})

test_that("plot draws the three types and returns x invisibly", {
  sel <- make_selection()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(sel))
  expect_no_error(plot(sel, type = "diagnostics"))
  expect_no_error(plot(sel, type = "pre", layout = "grid", seed = 2))
  out <- withVisible(plot(sel))
  expect_false(out$visible)
  expect_identical(out$value, sel)
  expect_error(plot(sel, type = "histogram"))
})

test_that("plot falls back to the post plot when diagnostics were not computed", {
  sel <- make_selection(diagnostics = FALSE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_message(plot(sel, type = "diagnostics"), "selected vs discarded instead")
  expect_no_message(plot(sel))
})

test_that("print shows a long call on one line without indentation", {
  sel <- make_selection()
  sel$call <- quote(rkn_select_criteria(g = g, criterion_primary = rkn_criterion_level("sex", "F"),
                                        criterion_secondary = rkn_criterion_cutoff("age", cutoff = 40, prefer = "below")))
  expect_output(print(sel), 'criterion_primary = rkn_criterion_level("sex", "F"), criterion_secondary', fixed = TRUE)
})
