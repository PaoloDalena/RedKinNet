adjacency <- function(g) as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE)) > 0

# Global selection of a network by proportion of the category, with the single
# individuals grouped as in the package.
global_selection <- function(g, cat, miss, k, p, ...) {
  comp <- igraph::components(g)
  single <- which(comp$csize[comp$membership] == 1)
  fams <- list()
  idx <- list()
  for (grp in .rkn_global_singletons(cat[single], miss[single], k)) {
    fams[[length(fams) + 1]] <- grp
    idx[[length(idx) + 1]] <- single[grp$members]
  }
  for (f in which(comp$csize > 1)) {
    v <- which(comp$membership == f)
    fams[[length(fams) + 1]] <- .rkn_global_options(adjacency(igraph::induced_subgraph(g, v)), cat[v], miss[v], k, ...)
    idx[[length(idx) + 1]] <- v
  }
  res <- .rkn_global_select(fams, k, proportion_distance(p))
  sets <- .rkn_global_realise(fams, res$choice)
  sel <- rep(FALSE, igraph::vcount(g))
  for (f in seq_along(sets)) sel[idx[[f]][sets[[f]]]] <- TRUE
  list(res = res, sel = sel)
}

random_network <- function(n_families) {
  pool <- c(atlas_connected(1, 5), atlas_connected(1, 1), atlas_connected(1, 1))
  repeat {
    parts <- pool[sample(length(pool), n_families)]
    for (i in seq_along(parts)) igraph::V(parts[[i]])$name <- paste0(i, "_", seq_len(igraph::vcount(parts[[i]])))
    g <- do.call(igraph::disjoint_union, unname(parts))
    if (igraph::vcount(g) <= 12) return(g)
  }
}

test_that("global selections reach the best proportion with the largest size, and count the selections that do", {
  set.seed(19)
  for (i in 1:30) {
    g <- random_network(sample(2:4, 1))
    n <- igraph::vcount(g)
    miss <- stats::rbinom(n, 1, 0.15)
    cat <- sample(0:1, n, replace = TRUE) * (1 - miss)
    k <- sample(0:2, 1)
    p <- stats::runif(1)
    best <- best_selections(g, cat, miss, k, p)
    out <- global_selection(g, cat, miss, k, p)
    expect_true(paste(which(out$sel), collapse = ",") %in% best)
    expect_equal(exp(out$res$log_count), length(best))
    expect_true(out$res$exact)
  }
})

test_that("global selections are drawn uniformly among the best selections", {
  set.seed(20)
  tested <- 0
  while (tested < 3) {
    g <- random_network(3)
    n <- igraph::vcount(g)
    miss <- stats::rbinom(n, 1, 0.15)
    cat <- sample(0:1, n, replace = TRUE) * (1 - miss)
    k <- sample(0:1, 1)
    p <- stats::runif(1)
    best <- best_selections(g, cat, miss, k, p)
    if (length(best) < 2 || length(best) > 8) next
    tested <- tested + 1
    drawn <- character(150 * length(best))
    for (d in seq_along(drawn)) drawn[d] <- paste(which(global_selection(g, cat, miss, k, p)$sel), collapse = ",")
    expect_setequal(unique(drawn), best)
    expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
  }
})

test_that("extreme targets are reached among more selections than double precision can count", {
  # 1100 related pairs of a woman and a man: 2^1100 maximum sets, one with all women,
  # choose(1100, 550) with half of them.
  pairs <- igraph::make_graph(seq_len(2200), directed = FALSE)
  cat <- rep(c(1, 0), 1100)
  all_women <- global_selection(pairs, cat, rep(0, 2200), 0, 1)
  expect_identical(sum(cat[all_women$sel]), 1100)
  expect_equal(all_women$res$log_count, 0)
  half <- global_selection(pairs, cat, rep(0, 2200), 0, 0.5)
  expect_identical(sum(cat[half$sel]), 550)
  expect_equal(half$res$log_count, lchoose(1100, 550))
})

test_that("families beyond the budget or too intricate are included with the options found, and flagged", {
  width_limit_only()
  set.seed(21)
  petersen <- igraph::make_graph("Petersen")
  cat <- rep(c(1, 0), 5)
  partial <- global_selection(petersen, cat, rep(0, 10), 0, 0.5, max_cells = 1)
  expect_false(partial$res$exact)
  expect_true(is.na(partial$res$log_count))
  expect_identical(sum(partial$sel), 4L)
  expect_identical(igraph::ecount(igraph::induced_subgraph(petersen, which(partial$sel))), 0)

  ring <- igraph::make_ring(7)
  fallback <- function() c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)
  fixed <- .rkn_global_options(adjacency(ring), rep(1, 7), rep(0, 7), 0, fallback = fallback, max_cells = 1, max_width = 1)
  expect_false(fixed$exact)
  expect_identical(fixed$options$c, 3)
  expect_error(.rkn_global_options(adjacency(ring), rep(1, 7), rep(0, 7), 0, max_cells = 1, max_width = 1),
               "too intricate")
})

test_that("a family gives up at most the size of its maximum sets, so a larger loss stays within the budget", {
  path <- igraph::make_ring(7, circular = FALSE)
  A <- adjacency(path)
  cat <- c(1, 0, 1, 1, 0, 0, 1)
  miss <- rep(0, 7)
  elim <- .rkn_elimination_order(A)
  limit <- .rkn_profile_budget(elim, cat, 4)[["cells"]]
  expect_gt(.rkn_profile_budget(elim, cat, 6)[["cells"]], limit)
  full <- .rkn_global_options(A, cat, miss, 6, max_cells = Inf)
  capped <- .rkn_global_options(A, cat, miss, 6, max_cells = limit)
  expect_true(capped$exact)
  key <- function(o) sort(paste(o$options$j, o$options$c, o$options$m, round(o$options$lc, 10)))
  expect_identical(key(capped), key(full))
})

test_that("a selection allowed to give up individuals reaches the composition of the trade-off table", {
  skip_on_cran()
  ids <- rkn_example$samples$id
  chk <- rkn_check_samples(rkn_example$kinship, "third", samples = ids)
  g <- rkn_graph_from_kinship(rkn_example$kinship, "third", samples = ids, exclude = chk$contamination$id,
                              covariates = rkn_example$samples, check = FALSE)
  tr <- rkn_tradeoff(g, "sex", c(F = 0.65), max_loss = 15)
  sel <- rkn_select_global(g, "sex", c(F = 0.65), max_loss = 15, seed = 1)
  expect_identical(attr(tr, "tradeoff")$n_inexact, 0L)
  expect_true(sel$global$exact)
  expect_equal(abs(sel$global$proportion - 0.65), min(tr$distance))
})

test_that("the memory of the global selection is checked before the tables are built", {
  b <- .rkn_global_singletons(rep(c(1, 0), 50), rep(0, 100), 2)
  expect_error(.rkn_global_select(b, 2, proportion_distance(0.5), max_cells = 10), "table cells")
})

test_that("the tables cover only the totals reachable with at most k individuals given up", {
  # Three individuals without relatives for each value from 1 to 50: giving up at most
  # three of them lowers the total by at most 150, whatever the values, while the
  # widths of the groups add up to 3825.
  fams <- list()
  for (v in 1:50) {
    t <- 0:3
    fams[[v]] <- list(alpha = 3, exact = TRUE, group = 3, members = 1:3,
                      options = list(j = t, c = (3 - t) * v, m = rep(0, 4), lc = lchoose(3, t), v = rep(NA_real_, 4)))
  }
  res <- .rkn_global_select(fams, 3, NULL, max_cells = 2000, draw = FALSE)
  total <- 3 * sum(1:50)
  expect_identical(dim(res$table), c(4L, 151L, 1L))
  expect_equal(res$c_base, total - 150)
  expect_equal(res$table[1, 151, 1], 0)
  # Giving up one individual of value v: three ways.
  expect_equal(res$table[2, total - (1:50) - res$c_base + 1, 1], rep(log(3), 50))
})

test_that("the work of the global tables is estimated before they are built", {
  b <- .rkn_global_singletons(rep(c(1, 0), 50), rep(0, 100), 2)
  expect_error(.rkn_global_select(b, 2, proportion_distance(0.5), max_work = 1, long = "stop"),
               class = "rkn_work_limit")
  expect_message(res <- .rkn_global_select(b, 2, proportion_distance(0.5), max_work = 1), "may take several minutes")
  expect_gt(res$work, 1)
  expect_silent(.rkn_global_select(b, 2, proportion_distance(0.5)))
})
