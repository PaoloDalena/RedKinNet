test_that("the three plotting functions run without error", {
  g <- make_network(31)
  res <- rkn_select_random(g, seed = 1, diagnostics = TRUE, verbose = FALSE)

  grDevices::pdf(NULL)
  expect_no_error(rkn_plot_pre(g))
  expect_no_error(rkn_plot_pre(g, layout = "gem", show_legend = TRUE))
  expect_no_error(rkn_plot_post(g, kept = res$kept))
  expect_no_error(rkn_plot_diagnostics(g, res = res))
  grDevices::dev.off()
})
