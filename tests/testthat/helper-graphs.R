# Shared fixtures and reference solvers for the test suite. The reference solvers
# are deliberately independent of igraph::largest_ivs(), which the package uses
# internally.

atlas_connected <- function(n_min = 1, n_max = 7) {
  out <- list()
  for (i in 1:1252) {
    g <- igraph::graph_from_atlas(i)
    n <- igraph::vcount(g)
    if (n < n_min || n > n_max || !igraph::is_connected(g)) next
    igraph::V(g)$name <- as.character(seq_len(n))
    out[[paste0("atlas_", i)]] <- g
  }
  out
}

# On CRAN (NOT_CRAN unset) exhaustive batteries run on a fixed one-in-ten subsample,
# so that every code path is still exercised within the CRAN time budget.
thin_on_cran <- function(x) {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) return(x)
  x[seq(1, length(x), by = 10)]
}

fam_df_from_graph <- function(g) {
  data.frame(
    .name = igraph::V(g)$name,
    name = igraph::V(g)$name,
    degree = as.integer(igraph::degree(g)),
    max_path_length = as.integer(igraph::eccentricity(g)),
    stringsAsFactors = FALSE
  )
}

add_covariates <- function(fam_df, seed) {
  set.seed(seed)
  n <- nrow(fam_df)
  fam_df$age <- sample(20:80, n, replace = TRUE)
  fam_df$sex <- sample(c("F", "M"), n, replace = TRUE)
  fam_df$grp <- sample(c("X", "Y", "Z"), n, replace = TRUE)
  fam_df
}

# All maximum independent sets by exhaustive enumeration (feasible up to ~15 nodes).
brute_force_mis <- function(g) {
  n <- igraph::vcount(g)
  ids <- igraph::V(g)$name
  adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  alpha <- 0
  sets <- list()
  for (k in 0:(2^n - 1)) {
    s <- which(bitwAnd(k, 2^(seq_len(n) - 1)) > 0)
    if (length(s) < alpha) next
    if (length(s) > 1 && any(adj[s, s] != 0)) next
    if (length(s) > alpha) {
      alpha <- length(s)
      sets <- list()
    }
    sets[[length(sets) + 1]] <- ids[s]
  }
  list(alpha = alpha, sets = sets)
}

# Independence number by integer linear programming (any graph size used in the tests).
ilp_alpha <- function(g) {
  n <- igraph::vcount(g)
  el <- igraph::as_edgelist(g, names = FALSE)
  if (nrow(el) == 0) return(n)
  a <- matrix(0, nrow(el), n)
  for (e in seq_len(nrow(el))) a[e, el[e, ]] <- 1
  sol <- lpSolve::lp("max", rep(1, n), a, rep("<=", nrow(el)), rep(1, nrow(el)), all.bin = TRUE)
  as.integer(round(sol$objval))
}

# Disjoint union of atlas families and random trees, with unique vertex names.
make_network <- function(seed, n_atlas = 15, tree_sizes = c(8, 12, 20)) {
  set.seed(seed)
  pool <- atlas_connected(1, 7)
  parts <- list()
  picks <- sample(seq_along(pool), n_atlas, replace = TRUE)
  for (j in seq_along(picks)) parts[[length(parts) + 1]] <- pool[[picks[j]]]
  for (s in tree_sizes) {
    parts[[length(parts) + 1]] <- igraph::sample_pa(s, m = 1, directed = FALSE)
  }
  g <- igraph::make_empty_graph(0, directed = FALSE)
  for (j in seq_along(parts)) {
    p <- parts[[j]]
    igraph::V(p)$name <- paste0("f", j, "_", seq_len(igraph::vcount(p)))
    g <- igraph::disjoint_union(g, p)
  }
  g
}

network_alpha <- function(g) {
  cc <- igraph::components(g)
  alpha <- 0
  for (k in seq_len(cc$no)) {
    sub <- igraph::induced_subgraph(g, which(cc$membership == k))
    alpha <- alpha + if (igraph::vcount(sub) <= 12) brute_force_mis(sub)$alpha else ilp_alpha(sub)
  }
  alpha
}

# Returns a character vector describing every violated invariant (empty if none).
selection_problems <- function(g, kept, discarded, alpha) {
  ids <- igraph::V(g)$name
  p <- character(0)
  if (anyDuplicated(kept) > 0) p <- c(p, "duplicated ids in kept")
  if (length(intersect(kept, discarded)) > 0) p <- c(p, "kept and discarded overlap")
  if (!setequal(union(kept, discarded), ids)) p <- c(p, "kept and discarded do not cover all nodes")
  if (length(kept) > 1 && igraph::ecount(igraph::induced_subgraph(g, kept)) > 0) {
    p <- c(p, "kept is not an independent set")
  }
  if (length(kept) != alpha) p <- c(p, sprintf("kept size %d differs from alpha %d", length(kept), alpha))
  p
}

make_selection <- function(diagnostics = TRUE) {
  g <- make_network(51)
  rkn_select_random(g, seed = 1, diagnostics = diagnostics)
}

# Synthetic kinship matrix of 200 samples: 30 parent-child pairs, a duplicate (S199,
# S200), a contaminated sample (S001) and 15 samples of a different ancestry (S181-S195).
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

# The network of the examples of the help pages, with the sex of each individual.
example_network <- function() {
  g <- igraph::graph_from_literal(P1 - C1 - P2, S1 - S2, T1 - T2 - T3 - T1,
                                  D1 - K1 - D2, D1 - K2 - D2, K1 - K2,
                                  D2 - K3 - D3, K1 - K3, K2 - K3, U1)
  sex <- c(P1 = "F", C1 = "M", P2 = "M", S1 = "M", S2 = "F", T1 = "M", T2 = "F", T3 = "M",
           D1 = "M", K1 = "F", D2 = "F", K2 = "M", K3 = "F", D3 = "M", U1 = "F")
  igraph::V(g)$sex <- unname(sex[igraph::V(g)$name])
  g
}

# Distance of a global selection from a target proportion p of the category, as in
# rkn_select_global().
proportion_distance <- function(p) {
  function(C, M, S) {
    D <- S - M
    d <- abs(C - p * D) / D
    d[D <= 0] <- Inf
    d
  }
}

# All the selections with the best proportion and, among them, the largest size.
best_selections <- function(g, cat, miss, k, p) {
  A <- .rkn_adjacency(g)
  n <- nrow(A)
  labels <- character(0)
  d <- numeric(0)
  size <- integer(0)
  for (s in 0:(2^n - 1)) {
    sel <- bitwAnd(s, 2^(seq_len(n) - 1)) > 0
    if (sum(sel) > 1 && any(A[sel, sel])) next
    labels <- c(labels, paste(which(sel), collapse = ","))
    size <- c(size, sum(sel))
    D <- sum(sel) - sum(miss[sel])
    d <- c(d, if (D > 0) abs(sum(cat[sel]) - p * D) / D else Inf)
  }
  keep <- size >= max(size) - k
  labels <- labels[keep]
  d <- d[keep]
  size <- size[keep]
  best <- d <= min(d) + 1e-12
  labels[best & size == max(size[best])]
}
