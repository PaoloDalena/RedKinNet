fam <- data.frame(
  .name = c("a", "b", "c", "d"),
  age = c(10, 20, 30, NA),
  sex = c("F", "M", "F", NA),
  grp = c("X", "X", "Y", "Z"),
  stringsAsFactors = FALSE
)

test_that("constructors return rkn_criterion objects and validate their inputs", {
  expect_s3_class(rkn_criterion_cutoff("age", 10), c("rkn_criterion_cutoff", "rkn_criterion"))
  expect_s3_class(rkn_criterion_level("sex", "F"), c("rkn_criterion_level", "rkn_criterion"))

  expect_error(rkn_criterion_cutoff(1, 10), "covariate")
  expect_error(rkn_criterion_cutoff("age", NA_real_), "cutoff")
  expect_error(rkn_criterion_cutoff("age", 10, prefer = "sideways"), "`prefer` must be \"above\" or \"below\"")
  expect_error(rkn_criterion_level("sex", c("F", "M")), "level")
})

test_that("cutoff scores the share of kept nodes on the preferred side", {
  expect_equal(rkn_criterion_cutoff("age", 15, "above")$score_fn(fam, c("a", "b", "c")), 2 / 3)
  expect_equal(rkn_criterion_cutoff("age", 15, "below")$score_fn(fam, c("a", "b", "c")), 1 / 3)
  expect_equal(rkn_criterion_cutoff("age", 15, "above")$score_fn(fam, "d"), 0)
})

test_that("level scores the share of kept nodes at the preferred level", {
  crit <- rkn_criterion_level("sex", "F")
  expect_equal(crit$score_fn(fam, c("a", "b", "c")), 2 / 3)
  expect_equal(crit$score_fn(fam, c("a", "d")), 1)
})

test_that("a missing covariate column raises an informative error", {
  expect_error(rkn_criterion_level("height", 1)$score_fn(fam, "a"), "not found")
})

test_that("print shows the criterion type and its parameters", {
  expect_output(print(rkn_criterion_level("sex", "F")), "<rkn_criterion: level>")
  expect_output(print(rkn_criterion_cutoff("age", 50, "below")), "prefer: below")
})
