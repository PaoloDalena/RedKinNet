fam <- data.frame(
  .name = c("a", "b", "c", "d"),
  age = c(10, 20, 30, NA),
  sex = c("F", "M", "F", NA),
  grp = c("X", "X", "Y", "Z"),
  stringsAsFactors = FALSE
)

test_that("constructors return rkn_criterion objects and validate their inputs", {
  expect_s3_class(rkn_criterion_target_mean("age", 10, 5), c("rkn_criterion_target_mean", "rkn_criterion"))
  expect_s3_class(rkn_criterion_cutoff("age", 10), c("rkn_criterion_cutoff", "rkn_criterion"))
  expect_s3_class(rkn_criterion_level("sex", "F"), c("rkn_criterion_level", "rkn_criterion"))
  expect_s3_class(rkn_criterion_category_balance("grp"), c("rkn_criterion_category_balance", "rkn_criterion"))

  expect_error(rkn_criterion_target_mean(1, 10, 5), "covariate")
  expect_error(rkn_criterion_target_mean("age", "10", 5), "target")
  expect_error(rkn_criterion_target_mean("age", 10, 0), "tolerance")
  expect_error(rkn_criterion_cutoff("age", NA_real_), "cutoff")
  expect_error(rkn_criterion_cutoff("age", 10, prefer = "sideways"))
  expect_error(rkn_criterion_level("sex", c("F", "M")), "level")
  expect_error(rkn_criterion_category_balance("grp", c(0.5, 0.5)), "named numeric")
  expect_error(rkn_criterion_category_balance("grp", c(X = 0.6, Y = 0.6)), "sum to 1")
})

test_that("target_mean decays linearly from the target and ignores missing values", {
  crit <- rkn_criterion_target_mean("age", target = 10, tolerance = 10)
  expect_equal(crit$score_fn(fam, c("a", "b")), 0.5)
  expect_equal(crit$score_fn(fam, "a"), 1)
  expect_equal(crit$score_fn(fam, "c"), 0)
  expect_equal(crit$score_fn(fam, c("a", "d")), 1)
  expect_equal(crit$score_fn(fam, character(0)), 0)
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

test_that("category_balance is one minus the total variation distance", {
  crit <- rkn_criterion_category_balance("grp", c(X = 0.5, Y = 0.5))
  expect_equal(crit$score_fn(fam, c("a", "b", "c")), 1 - 0.5 * (abs(2 / 3 - 0.5) + abs(1 / 3 - 0.5)))
  expect_equal(crit$score_fn(fam, c("a", "d")), 0.75)
  own_mix <- rkn_criterion_category_balance("grp")
  expect_equal(own_mix$score_fn(fam, c("a", "b", "c", "d")), 1)
})

test_that("a missing covariate column raises an informative error", {
  expect_error(rkn_criterion_level("height", 1)$score_fn(fam, "a"), "not found")
})

test_that("print shows the criterion type and its parameters", {
  expect_output(print(rkn_criterion_level("sex", "F")), "<rkn_criterion: level>")
  expect_output(print(rkn_criterion_category_balance("grp")), "target_props: NULL")
})
