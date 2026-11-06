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
    res <- rkn_select_random(g, seed = seed)
    expect_identical(selection_problems(g, res$kept, res$discarded, network_alpha(g)), character(0))
  }
})

test_that("network-level random selection is reproducible", {
  g <- make_network(11)
  expect_identical(rkn_select_random(g, seed = 3), rkn_select_random(g, seed = 3))
})

test_that("the selection is an rkn_selection object recording how it was made", {
  g <- make_network(13)
  sel <- rkn_select_random(g, seed = 1)
  expect_s3_class(sel, "rkn_selection")
  expect_identical(sel$method, "random")
  expect_null(sel$criteria)
  expect_identical(sel$seed, 1)
  expect_identical(sel$graph, g)
  expect_identical(sel$call, quote(rkn_select_random(g = g, seed = 1)))
  expect_identical(sel$family_sizes, as.integer(igraph::components(g)$csize))
  expect_true(sel$diagnostics)
})

test_that("diagnostics have the documented structure and agree with the selection", {
  g <- make_network(13)
  plain <- rkn_select_random(g, seed = 1, diagnostics = FALSE)
  expect_type(plain$kept, "character")
  expect_null(plain$node_diagnostics)
  expect_null(plain$family_summary)

  res <- rkn_select_random(g, seed = 1)
  expect_identical(res$kept, plain$kept)
  expect_identical(res$discarded, plain$discarded)
  expect_equal(nrow(res$node_diagnostics), igraph::vcount(g))
  expect_equal(nrow(res$family_summary), igraph::components(g)$no)
  nd <- res$node_diagnostics
  expect_setequal(nd$id[nd$is_selected], res$kept)
  expect_equal(sum(res$family_summary$n_kept), length(res$kept))
})

test_that("invalid inputs are rejected", {
  expect_error(rkn_select_random(list()), "igraph object")
  g_unnamed <- igraph::make_ring(5)
  expect_error(rkn_select_random(g_unnamed), "must be defined")
  expect_error(rkn_select_random(igraph::make_empty_graph(0)), "at least one vertex")
})

test_that("vertex attributes with reserved names are rejected instead of being mistaken for structural measures", {
  g <- make_network(7)
  igraph::V(g)$degree <- 1
  expect_error(rkn_select_random(g, seed = 1), 'rename "degree"')
  expect_error(rkn_select_criteria(g, rkn_criterion_level("degree", 1), seed = 1), 'rename "degree"')
  g <- make_network(7)
  igraph::V(g)$max_path_length <- 1
  expect_error(rkn_select_random(g, seed = 1), 'rename "max_path_length"')
})

test_that("covariates do not change a random selection", {
  for (s in 7:9) {
    g <- make_network(s)
    with_covariate <- g
    set.seed(s)
    igraph::V(with_covariate)$education <- sample(1:5, igraph::vcount(g), replace = TRUE)
    expect_identical(rkn_select_random(with_covariate, seed = 1)$kept, rkn_select_random(g, seed = 1)$kept)
  }
})

test_that("families left undecided by the structural rules are drawn uniformly by the exact engine", {
  petersen <- igraph::make_graph("Petersen")
  igraph::V(petersen)$name <- paste0("p", 1:10)
  ref <- brute_force_mis(petersen)
  drawn <- character(600)
  for (s in seq_along(drawn)) drawn[s] <- paste(sort(rkn_select_random(petersen, seed = s, diagnostics = FALSE)$kept), collapse = ",")
  expect_setequal(unique(drawn), vapply(ref$sets, function(x) paste(sort(x), collapse = ","), ""))
  expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
  expect_identical(rkn_select_random(petersen, seed = 1)$family_summary$mis_method, "exact")
})

test_that("large intricate families are resolved quickly with a maximum set", {
  grid <- igraph::make_lattice(c(6, 30))
  igraph::V(grid)$name <- paste0("g", seq_len(igraph::vcount(grid)))
  elapsed <- system.time(sel <- rkn_select_random(grid, seed = 1))[["elapsed"]]
  expect_lt(elapsed, 10)
  expect_length(sel$kept, 90)
  expect_identical(igraph::ecount(igraph::induced_subgraph(grid, sel$kept)), 0)
})

test_that("mis_method records how each family was resolved", {
  g <- igraph::disjoint_union(igraph::make_star(5, mode = "undirected"), igraph::make_full_graph(4),
                              igraph::make_graph("Petersen"))
  igraph::V(g)$name <- paste0("v", seq_len(igraph::vcount(g)))
  fs <- rkn_select_random(g, seed = 1)$family_summary
  expect_identical(fs$mis_method, c(NA, NA, "exact"))
})

test_that("individuals the search leaves without alternatives are decided without a choice", {
  # Family of six: the structural rules keep D3 and discard K3; the search resolves D1,
  # D2, K1 and K2, whose only maximum set is {D1, D2}. Family of five: two parents and
  # three children, resolved by enumeration, whose only maximum set is the children.
  g <- igraph::graph_from_literal(D1 - K1 - D2, D1 - K2 - D2, K1 - K2, D2 - K3 - D3, K1 - K3, K2 - K3,
                                  P1 - C1 - P2, P1 - C2 - P2, P1 - C3 - P2)
  igraph::V(g)$sex <- rep(c("M", "F"), length.out = 11)
  kept <- c("D1", "D2", "D3", "C1", "C2", "C3")
  discarded <- c("K1", "K2", "K3", "P1", "P2")
  selectors <- list(function(d) rkn_select_random(g, seed = 1, diagnostics = d),
                    function(d) rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1, diagnostics = d))
  for (select in selectors) {
    sel <- select(TRUE)
    expect_identical(sel$kept, select(FALSE)$kept)
    nd <- sel$node_diagnostics
    expect_true(all(nd$pre_class[nd$id %in% kept] == "kept_apriori"))
    expect_true(all(nd$pre_class[nd$id %in% discarded] == "discarded_apriori"))
    expect_setequal(sel$family_summary$mis_method, c("exact", "enumeration"))
    expect_identical(sel$family_summary$n_candidate_complex, c(0L, 0L))
    expect_setequal(sel$kept, kept)
  }
})

test_that("uniform selections draw every maximum set of a family with the same probability", {
  path6 <- igraph::make_ring(6, circular = FALSE)
  igraph::V(path6)$name <- paste0("p", 1:6)
  drawn <- character(800)
  for (s in seq_along(drawn)) {
    drawn[s] <- paste(sort(rkn_select_random(path6, seed = s, diagnostics = FALSE, uniform = TRUE)$kept), collapse = ",")
  }
  expect_setequal(unique(drawn), c("p1,p3,p5", "p1,p3,p6", "p1,p4,p6", "p2,p4,p6"))
  expect_gt(stats::chisq.test(table(drawn))$p.value, 0.001)
})

test_that("uniform diagnostics give exact inclusion probabilities and decision types", {
  g <- igraph::disjoint_union(igraph::make_ring(6, circular = FALSE), igraph::make_graph("Petersen"),
                              igraph::make_star(5, mode = "undirected"))
  ids <- paste0("v", seq_len(igraph::vcount(g)))
  igraph::V(g)$name <- ids
  sel <- rkn_select_random(g, seed = 1, uniform = TRUE)
  expect_true(sel$uniform)
  expect_true(all(sel$family_summary$uniform))
  expect_identical(sel$family_summary$selector, rep("uniform", 3))
  nd <- sel$node_diagnostics
  comp <- igraph::components(g)$membership
  for (f in 1:3) {
    sub <- igraph::induced_subgraph(g, which(comp == f))
    sets <- brute_force_mis(sub)$sets
    members <- igraph::V(sub)$name
    expected <- numeric(length(members))
    for (s in sets) expected <- expected + (members %in% s)
    expected <- expected / length(sets)
    rows <- match(members, nd$id)
    expect_equal(nd$inclusion_probability[rows], expected)
    expect_identical(nd$pre_class[rows],
                     ifelse(expected == 1, "kept_apriori", ifelse(expected == 0, "discarded_apriori", "candidate_complex")))
    expect_equal(sel$family_summary$n_max_sets[nd$family_id[rows[1]]], length(sets))
  }
  expect_output(print(sel), "Uniform random selection")
  expect_output(print(sel), "Families sampled uniformly among their maximum sets: 3 of 3")
  expect_output(print(sel), "Maximum sets: 20 in total; families with more than one: 2 of 3", fixed = TRUE)
  expect_output(print(summary(sel)), "in every maximum set")
  expect_true(all(is.na(rkn_select_random(g, seed = 1)$family_summary$n_max_sets)))

  pairs <- igraph::make_graph(seq_len(2 * 30), directed = FALSE)
  igraph::V(pairs)$name <- paste0("q", seq_len(60))
  expect_output(print(rkn_select_random(pairs, seed = 1, uniform = TRUE)),
                "Maximum sets: 1.1e+09 in total; families with more than one: 30 of 30", fixed = TRUE)
  # Each member of a pair is selected with probability 1/2: 60 x 1/4 = 15 differences.
  expect_output(print(rkn_select_random(pairs, seed = 1, uniform = TRUE)),
                "Two uniform draws differ on average in 15.0 of the 30 selected individuals (50.0%)", fixed = TRUE)
  expect_false(any(grepl("Two uniform draws", capture.output(print(rkn_select_random(pairs, seed = 1))))))
})

test_that("families too intricate for uniform sampling are selected by the rules and reported", {
  petersen <- igraph::make_graph("Petersen")
  igraph::V(petersen)$name <- paste0("p", 1:10)
  fam_df <- fam_df_from_graph(petersen)
  out <- .rkn_select_uniform_family(fam_df, petersen, seed = 1, return_diagnostics = TRUE, max_width = 1)
  expect_false(out$uniform)
  expect_length(out$kept, 4)
  expect_true(all(is.na(out$node_diagnostics$inclusion_probability)))
  expect_false(rkn_select_random(petersen, seed = 1)$family_summary$uniform)
})
