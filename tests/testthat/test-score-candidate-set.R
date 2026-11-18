fam <- data.frame(
  .name = c("a", "b", "c"),
  age = c(30, 50, 70),
  sex = c("F", "F", "M"),
  stringsAsFactors = FALSE
)

test_that(".rkn_score_candidate_set returns one row with primary and secondary scores", {
  out <- .rkn_score_candidate_set(fam, c("a", "c"), rkn_criterion_level("sex", "F"))
  expect_equal(nrow(out), 1)
  expect_named(out, c("primary_type", "primary_score", "secondary_type", "secondary_score"))
  expect_equal(out$primary_score, 0.5)
  expect_true(is.na(out$secondary_score))

  out2 <- .rkn_score_candidate_set(fam, "b", rkn_criterion_level("sex", "F"),
                                  rkn_criterion_cutoff("age", 40, "above"))
  expect_equal(out2$secondary_score, 1)
  expect_identical(out2$secondary_type, "cutoff")
})

test_that(".rkn_score_candidate_set validates its criteria", {
  expect_error(.rkn_score_candidate_set(fam, "a", "x"), "criterion_primary")
  expect_error(.rkn_score_candidate_set(fam, "a", rkn_criterion_level("sex", "F"), "x"), "criterion_secondary")
})

test_that("candidate choice follows the primary score, then the secondary, then a seeded draw", {
  select <- .rkn_select_among_candidates
  level_f <- rkn_criterion_level("sex", "F")
  over_40 <- rkn_criterion_cutoff("age", 40, "above")

  expect_identical(select(list("c", "a"), fam, criterion_primary = level_f)$chosen_set, "a")
  expect_identical(select(list("a", "b"), fam, criterion_primary = level_f,
                          criterion_secondary = over_40)$chosen_set, "b")

  tied <- list("a", "b")
  first <- select(tied, fam, criterion_primary = level_f, seed = 5)$chosen_set
  for (i in 1:5) {
    expect_identical(select(tied, fam, criterion_primary = level_f, seed = 5)$chosen_set, first)
  }
  expect_identical(select(list("c"), fam, criterion_primary = level_f)$chosen_set, "c")
  expect_error(select(list(), fam, criterion_primary = level_f), "non-empty")
})

test_that("fixed kept nodes enter the score of every candidate", {
  select <- .rkn_select_among_candidates
  fam_na <- data.frame(.name = c("a", "b", "c"), sex = c(NA, "F", "M"), stringsAsFactors = FALSE)
  level_f <- rkn_criterion_level("sex", "F")
  # Alone, a (missing value) and c (a man) both score 0; with b, a scores 1 and c 1/2.
  for (s in 1:5) {
    expect_identical(select(list("a", "c"), fam_na, fixed_kept_ids = "b", criterion_primary = level_f,
                            seed = s)$chosen_set, "a")
  }
})
