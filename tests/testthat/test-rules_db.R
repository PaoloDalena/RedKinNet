test_that("every connected graph with at most 5 nodes matches exactly one rule", {
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

test_that("every rule in rules_db is reachable by at least one connected graph", {
  graphs <- atlas_connected(1, 5)
  used <- character(0)
  for (nm in names(graphs)) {
    rule <- .rkn_find_matching_rule(fam_df_from_graph(graphs[[nm]]), rules_db)
    used <- c(used, rule$name)
  }
  all_names <- character(0)
  for (rule in rules_db) all_names <- c(all_names, rule$name)
  expect_identical(setdiff(all_names, used), character(0))
})
