kinship_example <- function() {
  K <- matrix(-0.01, 6, 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
  set_pair <- function(a, b, phi) {
    K[a, b] <<- phi
    K[b, a] <<- phi
  }
  set_pair("A", "B", 0.25)
  set_pair("C", "D", 0.06)
  set_pair("A", "C", 0.0442)
  set_pair("E", "F", 0.12)
  diag(K) <- 0.5
  K
}

edge_set <- function(g) {
  e <- igraph::as_edgelist(g)
  sort(paste(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2])))
}

all_pairs <- function(K) {
  ij <- which(upper.tri(K), arr.ind = TRUE)
  data.frame(ID1 = rownames(K)[ij[, 1]], ID2 = rownames(K)[ij[, 2]], Kinship = K[ij])
}

test_that("a matrix gives one vertex per sample and an edge for each pair at or above the threshold", {
  K <- kinship_example()
  g <- rkn_graph_from_kinship(K, threshold = 0.0442)
  expect_false(igraph::is_directed(g))
  expect_identical(igraph::V(g)$name, LETTERS[1:6])
  expect_identical(edge_set(g), c("A B", "A C", "C D", "E F"))
  expect_equal(sort(igraph::E(g)$kinship), sort(c(0.25, 0.0442, 0.06, 0.12)))
  expect_identical(g$kinship_threshold, 0.0442)
  expect_identical(edge_set(rkn_graph_from_kinship(K, threshold = 0.1)), c("A B", "E F"))
})

test_that("threshold labels select the KING cutoffs", {
  K <- kinship_example()
  expect_identical(rkn_graph_from_kinship(K, "third")$kinship_threshold, 0.0442)
  expect_identical(rkn_graph_from_kinship(K, "second")$kinship_threshold, 0.0884)
  expect_identical(rkn_graph_from_kinship(K, "first")$kinship_threshold, 0.177)
  expect_identical(rkn_graph_from_kinship(K, "duplicate")$kinship_threshold, 0.354)
  expect_identical(edge_set(rkn_graph_from_kinship(K, "second")), c("A B", "E F"))
  expect_identical(edge_set(rkn_graph_from_kinship(K, "first")), "A B")
  expect_error(rkn_graph_from_kinship(K, "cousins"))
  expect_error(rkn_graph_from_kinship(K, c(0.1, 0.2)), "threshold must be")
})

test_that("a complete pair table gives the same network as the matrix", {
  K <- kinship_example()
  from_matrix <- rkn_graph_from_kinship(K, "third")
  from_pairs <- expect_no_message(rkn_graph_from_kinship(all_pairs(K), "third"))
  expect_identical(sort(igraph::V(from_pairs)$name), igraph::V(from_matrix)$name)
  expect_identical(edge_set(from_pairs), edge_set(from_matrix))
})

test_that("pair tables are read with KING and PLINK 2 column names", {
  pairs <- all_pairs(kinship_example())
  plink <- data.frame(`#FID1` = 0, IID1 = pairs$ID1, FID2 = 0, IID2 = pairs$ID2, KINSHIP = pairs$Kinship,
                      check.names = FALSE)
  expect_identical(edge_set(rkn_graph_from_kinship(plink, "third")), edge_set(rkn_graph_from_kinship(pairs, "third")))
  names(pairs) <- c("a", "b", "phi")
  expect_error(rkn_graph_from_kinship(pairs, "third"), "two identifier columns")
})

test_that("samples restores the unrelated individuals missing from a pair table", {
  related <- data.frame(ID1 = c("A", "C"), ID2 = c("B", "D"), Kinship = c(0.25, 0.06))
  expect_message(g <- rkn_graph_from_kinship(related, "third"), "samples without relatives")
  expect_identical(igraph::V(g)$name, c("A", "C", "B", "D"))
  g <- expect_no_message(rkn_graph_from_kinship(related, "third", samples = LETTERS[1:6]))
  expect_identical(igraph::V(g)$name, LETTERS[1:6])
  expect_identical(igraph::degree(g, c("E", "F")), c(E = 0, F = 0))
  expect_error(rkn_graph_from_kinship(related, "third", samples = c("A", "B", "C")), "missing from `samples`")
})

test_that("exclude removes samples and their relationships", {
  g <- rkn_graph_from_kinship(kinship_example(), "third", exclude = "A")
  expect_identical(igraph::V(g)$name, LETTERS[2:6])
  expect_identical(edge_set(g), c("C D", "E F"))
  expect_error(rkn_graph_from_kinship(kinship_example(), "third", exclude = "Z"), "absent from the data")
})

test_that("covariates become vertex attributes, with NA for missing samples", {
  covariates <- data.frame(sample = c("F", "A", "B"), sex = c("F", "M", "F"), age = c(30, 50, 20))
  expect_message(g <- rkn_graph_from_kinship(kinship_example(), "third", covariates = covariates,
                                             covariates_id = "sample"),
                 "missing for 3 of 6")
  expect_identical(igraph::V(g)$sex, c("M", "F", NA, NA, NA, "F"))
  expect_identical(igraph::V(g)$age, c(50, 20, NA, NA, NA, 30))
  expect_error(rkn_graph_from_kinship(kinship_example(), "third", covariates = covariates), "named id")
  covariates$sample[2] <- "F"
  expect_error(rkn_graph_from_kinship(kinship_example(), "third", covariates = covariates, covariates_id = "sample"),
               "duplicated")
})

test_that("invalid matrices and missing estimates are handled", {
  K <- kinship_example()
  expect_error(rkn_graph_from_kinship(unname(K), "third"), "row and column names")
  asym <- K
  asym["A", "B"] <- 0.1
  expect_error(rkn_graph_from_kinship(asym, "third"), "symmetric")
  K["E", "F"] <- K["F", "E"] <- NA
  expect_warning(g <- rkn_graph_from_kinship(K, "third"), "1 pairs have a missing kinship")
  expect_identical(edge_set(g), c("A B", "A C", "C D"))
  expect_error(rkn_graph_from_kinship(list(), "third"), "kinship matrix or a data frame")
})

test_that("the network can be passed to the selectors", {
  g <- rkn_graph_from_kinship(kinship_example(), "third")
  sel <- rkn_select_random(g, seed = 1)
  expect_s3_class(sel, "rkn_selection")
  expect_identical(sort(c(sel$kept, sel$discarded)), LETTERS[1:6])
})
