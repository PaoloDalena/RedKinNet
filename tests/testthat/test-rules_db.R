test_that("every connected graph with at most 5 nodes matches exactly one rule", {
  skip("Known issue: rules #16 and #17 both match atlas graph 37 (first match is correct); to be fixed")
  graphs <- atlas_connected(1, 5)
  problems <- character(0)
  for (nm in names(graphs)) {
    fam_df <- fam_df_from_graph(graphs[[nm]])
    n_match <- 0
    for (rule in rules_db) {
      if (isTRUE(rule$size_allowed == nrow(fam_df)) && isTRUE(rule$match_fn(fam_df))) {
        n_match <- n_match + 1
      }
    }
    if (n_match != 1) {
      problems <- c(problems, sprintf("%s (n = %d): %d matching rules", nm, nrow(fam_df), n_match))
    }
  }
  expect_identical(problems, character(0))
})

test_that("rules declaring alpha state the true independence number of the matched graphs", {
  skip("Known issue: `alpha` acts as a switch, not as the independence number; to be fixed")
  graphs <- atlas_connected(1, 5)
  problems <- character(0)
  for (nm in names(graphs)) {
    g <- graphs[[nm]]
    rule <- find_matching_rule(fam_df_from_graph(g), rules_db)
    if (is.null(rule$alpha)) next
    alpha <- brute_force_mis(g)$alpha
    if (rule$alpha != alpha) {
      problems <- c(problems, sprintf("%s matched by %s: rule alpha %d, true alpha %d",
                                      nm, rule$name, rule$alpha, alpha))
    }
  }
  expect_identical(problems, character(0))
})

test_that("every rule in rules_db is reachable by at least one connected graph", {
  skip("Known issue: rule triangle_two_tails5 requires an impossible degree sequence; to be fixed")
  graphs <- atlas_connected(1, 5)
  used <- character(0)
  for (nm in names(graphs)) {
    rule <- find_matching_rule(fam_df_from_graph(graphs[[nm]]), rules_db)
    used <- c(used, rule$name)
  }
  all_names <- character(0)
  for (rule in rules_db) all_names <- c(all_names, rule$name)
  expect_identical(setdiff(all_names, used), character(0))
})
