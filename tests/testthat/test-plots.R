test_that("rkn_plot_pre draws with each layout, colouring and title option", {
  g <- make_network(31)

  grDevices::pdf(NULL)
  expect_no_error(rkn_plot_pre(g))
  expect_no_error(rkn_plot_pre(g, layout = "grid", seed = 2, colour_by = "family", main = ""))
  expect_no_error(rkn_plot_pre(g, vertex_size = 2, edge_width = 2, vertex_label = TRUE,
                               show_legend = FALSE, main = "A custom title"))
  grDevices::dev.off()
})

test_that("plot draws every type with both palettes", {
  sel <- make_selection()

  grDevices::pdf(NULL)
  for (type in c("post", "diagnostics", "pre")) {
    for (palette in c("colorblind", "classic")) {
      expect_no_error(plot(sel, type = type, palette = palette))
    }
  }
  grDevices::dev.off()
})

test_that("the plots restore only the margins, so panels can be combined", {
  sel <- make_selection()

  grDevices::pdf(NULL)
  graphics::par(mfrow = c(1, 3))
  mar <- graphics::par("mar")
  plot(sel, type = "pre")
  expect_equal(graphics::par("mfg"), c(1, 1, 1, 3))
  plot(sel)
  expect_equal(graphics::par("mfg"), c(1, 2, 1, 3))
  plot(sel, type = "diagnostics")
  expect_equal(graphics::par("mfg"), c(1, 3, 1, 3))
  expect_equal(graphics::par("mar"), mar)
  grDevices::dev.off()
})

test_that("the post plot colours selected and discarded individuals with the palette", {
  sel <- make_selection()
  kept <- igraph::V(sel$graph)$name %in% sel$kept

  grDevices::pdf(NULL)
  cb <- .rkn_plot_post(sel$graph, kept = sel$kept)
  classic <- .rkn_plot_post(sel$graph, kept = sel$kept, palette = "classic")
  grDevices::dev.off()
  expect_true(all(cb$fill[kept] == "#009E73"))
  expect_true(all(cb$fill[!kept] == "#999999"))
  expect_true(all(classic$fill[kept] == "#1B9E77"))
  expect_true(all(classic$fill[!kept] == "#B3B3B3"))
})

test_that("the diagnostics plot fills selected individuals and leaves discarded ones hollow", {
  sel <- make_selection()
  kept <- igraph::V(sel$graph)$name %in% sel$kept

  grDevices::pdf(NULL)
  outcome <- .rkn_plot_diagnostics(sel$graph, res = sel)
  plain <- .rkn_plot_diagnostics(sel$graph, res = sel, show_outcome = FALSE)
  grDevices::dev.off()
  expect_true(all(outcome$fill[!kept] == "white"))
  expect_equal(outcome$frame_col[!kept], plain$fill[!kept])
  expect_equal(outcome$fill[kept], plain$fill[kept])
  expect_false(any(plain$fill == "white"))
})

test_that("min_family_size draws only the families with at least that many members", {
  sel <- make_selection()
  sizes <- igraph::components(sel$graph)$csize
  n_large <- sum(sizes[sizes >= 8])
  expect_lt(n_large, igraph::vcount(sel$graph))

  grDevices::pdf(NULL)
  post <- .rkn_plot_post(sel$graph, kept = sel$kept, min_family_size = 8)
  diagnostics <- .rkn_plot_diagnostics(sel$graph, res = sel, min_family_size = 8)
  expect_no_error(plot(sel, type = "pre", min_family_size = 8))
  expect_error(plot(sel, min_family_size = 100), "No family has at least 100 members")
  grDevices::dev.off()
  expect_length(post$fill, n_large)
  expect_length(diagnostics$fill, n_large)
})
