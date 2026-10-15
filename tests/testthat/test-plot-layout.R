test_that("the layout is reproducible for a given seed and changes with the seed", {
  g <- make_network(31)
  for (layout in c("components", "grid")) {
    xy <- .rkn_layout(g, layout, seed = 1)
    expect_identical(.rkn_layout(g, layout, seed = 1), xy)
    expect_false(isTRUE(all.equal(.rkn_layout(g, layout, seed = 2), xy)))
    expect_equal(max(abs(xy)), 1)
  }
})

test_that("auto uses the components layout up to 500 individuals and the grid layout above", {
  pairs <- function(n_pairs) igraph::make_graph(seq_len(2 * n_pairs), directed = FALSE)
  g500 <- pairs(250)
  g502 <- pairs(251)
  expect_identical(.rkn_layout(g500, "auto", seed = 1), .rkn_layout(g500, "components", seed = 1))
  expect_identical(.rkn_layout(g502, "auto", seed = 1), .rkn_layout(g502, "grid", seed = 1))
  expect_false(identical(.rkn_layout(g500, "components", seed = 1), .rkn_layout(g500, "grid", seed = 1)))
})

test_that("the grid layout gives each family a region that no other family enters", {
  g <- make_network(31)
  xy <- .rkn_layout(g, "grid", seed = 1)
  fams <- split(seq_len(igraph::vcount(g)), igraph::components(g)$membership)
  box <- matrix(0, length(fams), 4)
  for (k in seq_along(fams)) box[k, ] <- c(range(xy[fams[[k]], 1]), range(xy[fams[[k]], 2]))
  overlap <- FALSE
  for (i in seq_along(fams)) {
    for (j in seq_along(fams)) {
      if (i < j) {
        overlap <- overlap || (box[i, 1] <= box[j, 2] && box[j, 1] <= box[i, 2] &&
                                 box[i, 3] <= box[j, 4] && box[j, 3] <= box[i, 4])
      }
    }
  }
  expect_false(overlap)
})
