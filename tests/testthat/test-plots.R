test_that("the three plotting functions run without error", {
  g <- make_network(31)
  res <- rkn_select_random(g, seed = 1, diagnostics = TRUE, verbose = FALSE)

  grDevices::pdf(NULL)
  expect_no_error(rkn_plot_pre(g))
  expect_no_error(rkn_plot_pre(g, layout = "grid", seed = 2, show_legend = TRUE))
  expect_no_error(rkn_plot_post(g, kept = res$kept))
  expect_no_error(rkn_plot_diagnostics(g, res = res))
  grDevices::dev.off()
})

test_that("the plotting functions restore only the margins, so panels can be combined", {
  g <- make_network(31)
  res <- rkn_select_random(g, seed = 1, diagnostics = TRUE, verbose = FALSE)

  grDevices::pdf(NULL)
  graphics::par(mfrow = c(1, 3))
  mar <- graphics::par("mar")
  rkn_plot_pre(g)
  expect_equal(graphics::par("mfg"), c(1, 1, 1, 3))
  rkn_plot_post(g, kept = res$kept)
  expect_equal(graphics::par("mfg"), c(1, 2, 1, 3))
  rkn_plot_diagnostics(g, res = res)
  expect_equal(graphics::par("mfg"), c(1, 3, 1, 3))
  expect_equal(graphics::par("mar"), mar)
  grDevices::dev.off()
})
