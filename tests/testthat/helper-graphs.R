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
