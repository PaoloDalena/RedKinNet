synthetic_cohort <- function() {
  set.seed(1)
  ids <- sprintf("S%03d", 1:200)
  K <- matrix(rnorm(200^2, -0.015, 0.004), 200, 200, dimnames = list(ids, ids))
  K[upper.tri(K)] <- t(K)[upper.tri(K)]
  ancestry <- 181:195
  K[ancestry, -ancestry] <- rnorm(15 * 185, -0.15, 0.03)
  K[-ancestry, ancestry] <- t(K[ancestry, -ancestry])
  for (i in seq(2, 60, by = 2)) K[i, i + 1] <- K[i + 1, i] <- 0.25
  K["S199", "S200"] <- K["S200", "S199"] <- 0.49
  K["S001", -1] <- K[-1, "S001"] <- rnorm(199, 0.045, 0.01)
  diag(K) <- 0.5
  K
}

pairs_of <- function(K) {
  ij <- which(upper.tri(K), arr.ind = TRUE)
  data.frame(ID1 = rownames(K)[ij[, 1]], ID2 = rownames(K)[ij[, 2]], Kinship = K[ij])
}

test_that("the check flags duplicates, contamination and different ancestry in a synthetic cohort", {
  chk <- rkn_check_samples(synthetic_cohort(), threshold = "third")
  expect_s3_class(chk, "rkn_sample_check")
  expect_identical(chk$duplicates[, c("id1", "id2")], data.frame(id1 = "S199", id2 = "S200"))
  expect_identical(chk$contamination$id, "S001")
  expect_identical(chk$contamination$n_relatives, 105L)
  expect_identical(chk$largest_component, 121L)
  expect_identical(chk$contamination$largest_component_without, 2L)
  ancestry <- chk$samples$id[chk$samples$flag == "possible different ancestry"]
  expect_identical(ancestry, sprintf("S%03d", 181:195))
  expect_true(chk$mean_kinship_available)
})

test_that("relatives are counted in total and by degree", {
  ids <- LETTERS[1:5]
  K <- matrix(-0.01, 5, 5, dimnames = list(ids, ids))
  K["A", "B"] <- K["B", "A"] <- 0.5
  K["A", "C"] <- K["C", "A"] <- 0.25
  K["A", "D"] <- K["D", "A"] <- 0.1
  K["A", "E"] <- K["E", "A"] <- 0.0442
  diag(K) <- 0.5
  s <- rkn_check_samples(K, "second")$samples
  expect_identical(s$n_relatives, c(3L, 1L, 1L, 1L, 0L))
  expect_identical(s$n_duplicate, c(1L, 1L, 0L, 0L, 0L))
  expect_identical(s$n_first, c(1L, 0L, 1L, 0L, 0L))
  expect_identical(s$n_second, c(1L, 0L, 0L, 1L, 0L))
  expect_identical(s$n_third, c(1L, 0L, 0L, 0L, 1L))
  # The pair with A reaches the third-degree cutoff, so it is left out of the mean.
  expect_equal(s$mean_kinship[5], -0.01)
})

test_that("a complete pair table gives the same check as the matrix", {
  K <- synthetic_cohort()
  from_matrix <- rkn_check_samples(K, "third")
  from_pairs <- expect_no_message(rkn_check_samples(pairs_of(K), "third", samples = rownames(K)))
  expect_equal(from_pairs$samples, from_matrix$samples)
  expect_equal(from_pairs$contamination, from_matrix$contamination)
})

test_that("with related pairs only, the mean kinship is unavailable and contamination uses the relatives", {
  K <- synthetic_cohort()
  related <- pairs_of(K)
  related <- related[related$Kinship >= 0.0442, ]
  expect_message(rkn_check_samples(related, "third"), "samples without relatives")
  chk <- rkn_check_samples(related, "third", samples = rownames(K))
  expect_false(chk$mean_kinship_available)
  expect_true(all(is.na(chk$samples$mean_kinship)))
  expect_identical(chk$contamination$id, "S001")
  expect_false(any(chk$samples$flag == "possible different ancestry"))
  unlisted <- setdiff(rownames(K), c(related$ID1, related$ID2))
  expect_true(length(unlisted) > 0)
  expect_true(all(chk$samples$n_relatives[chk$samples$id %in% unlisted] == 0L))
  expect_output(print(chk), "Ancestry not assessed")
})

test_that("robust z-scores fall back on the mean absolute deviation when the MAD is zero", {
  expect_equal(.rkn_robust_z(c(1, 2, 3, 4, 5)), (c(1, 2, 3, 4, 5) - 3) / stats::mad(1:5))
  expect_equal(.rkn_robust_z(c(0, 0, 0, 0, 5))[5], 5 / (sqrt(pi / 2) * 1))
  expect_identical(.rkn_robust_z(c(2, 2, 2)), c(0, 0, 0))
})

test_that("the print shows each section and the flagged samples can be excluded", {
  K <- synthetic_cohort()
  chk <- rkn_check_samples(K, "third")
  expect_output(print(chk), "Possible duplicates or monozygotic twins")
  expect_output(print(chk), "S001")
  expect_output(print(chk), "Possible different ancestry .*: 15 sample")
  capture.output(out <- withVisible(print(chk)))
  expect_false(out$visible)
  g <- rkn_graph_from_kinship(K, "third", exclude = chk$contamination$id)
  expect_identical(max(igraph::components(g)$csize), 2)
})

test_that("true relatives do not make the members of a small cohort look contaminated", {
  set.seed(2)
  ids <- sprintf("S%02d", 1:40)
  K <- matrix(rnorm(40^2, -0.015, 0.004), 40, 40, dimnames = list(ids, ids))
  K[upper.tri(K)] <- t(K)[upper.tri(K)]
  for (i in c(1, 3, 5)) K[i, i + 1] <- K[i + 1, i] <- 0.25
  diag(K) <- 0.5
  chk <- rkn_check_samples(K, "third")
  expect_identical(nrow(chk$contamination), 0L)
  expect_true(all(chk$samples$flag == ""))
})
