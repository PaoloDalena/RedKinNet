create_family_graph <- function(n, type = c("path", "cycle", "star", "clique", "random", "two_cycles")) {
  type <- match.arg(type)

  if (type == "path") {
    edge_pairs <- cbind(1:(n - 1), 2:n)
    g <- igraph::make_graph(edges = as.vector(t(edge_pairs)), n = n, directed = FALSE)

  } else if (type == "cycle") {
    g <- igraph::make_ring(n, directed = FALSE)

  } else if (type == "star") {
    g <- igraph::make_star(n, mode = "undirected", center = 1)

  } else if (type == "clique") {
    g <- igraph::make_full_graph(n, directed = FALSE)

  } else if (type == "two_cycles") {
    # figura a "8": due cicli che condividono il vertice 1
    if (n < 7) {
      g <- igraph::make_ring(n, directed = FALSE)
    } else {
      k <- 4  # primo ciclo: 1..4
      # ciclo1 su 1..k
      cyc1_pairs <- cbind(1:k, c(2:k, 1))
      e1 <- as.vector(t(cyc1_pairs))

      # ciclo2 su 1 e (k+1..n)
      nodes2 <- c(1, (k + 1):n)
      cyc2_pairs <- cbind(nodes2, c(nodes2[-1], nodes2[1]))
      e2 <- as.vector(t(cyc2_pairs))

      g <- igraph::make_graph(edges = c(e1, e2), n = n, directed = FALSE)
      g <- igraph::simplify(g)
    }

  } else if (type == "random") {
    p <- min(0.5, max(0.2, 2 / n))
    g <- igraph::sample_gnp(n, p, directed = FALSE, loops = FALSE)
    g <- igraph::simplify(g)

    # se non connesso, collega componenti in modo sicuro
    if (!igraph::is_connected(g)) {
      memb <- igraph::components(g)$membership
      # rappresentante per componente (sempre valido)
      reps <- unname(tapply(seq_len(n), memb, function(idx) idx[1]))
      reps <- as.integer(reps)

      if (length(reps) > 1) {
        for (i in seq_len(length(reps) - 1)) {
          g <- igraph::add_edges(g, c(reps[i], reps[i + 1]))
        }
      }
    }

    g <- igraph::simplify(g)
  }

  g
}

make_synthetic_network <- function(n6 = 50, n7 = 20, n8 = 10, n9 = 5, n10 = 3, seed = 123) {
  set.seed(seed)

  type_pool <- c("path", "cycle", "star", "clique", "random", "two_cycles")
  type_prob <- c(0.01, 0.01, 0.01, 0.01, 0.95, 0.01)

  make_block <- function(size, k) {
    if (k <= 0) {
      return(list(graphs = list(), meta = list()))
    }
    family_list <- vector("list", k)
    meta_list <- vector("list", k)

    for (i in seq_len(k)) {
      type <- sample(type_pool, 1, prob = type_prob)
      g <- create_family_graph(size, type = type)
      family_list[[i]] <- g
      meta_list[[i]] <- list(size = size, type = type)
    }
    list(graphs = family_list, meta = meta_list)
  }

  blocks <- list(
    make_block(6,  n6),
    make_block(7,  n7),
    make_block(8,  n8),
    make_block(9,  n9),
    make_block(10, n10)
  )

  graphs <- unlist(lapply(blocks, `[[`, "graphs"), recursive = FALSE)
  meta   <- unlist(lapply(blocks, `[[`, "meta"),   recursive = FALSE)

  if (length(graphs) == 0) stop("Hai richiesto 0 famiglie totali (tutti i n* sono 0).")

  g_all <- do.call(igraph::disjoint_union, graphs)
  g_all <- igraph::simplify(g_all)

  family_sizes <- purrr::map_int(meta, "size")
  family_ids <- rep(seq_along(family_sizes), times = family_sizes)

  igraph::V(g_all)$family_id <- family_ids
  igraph::V(g_all)$local_id  <- ave(family_ids, family_ids, FUN = seq_along)
  igraph::V(g_all)$name <- paste0("F", igraph::V(g_all)$family_id, "_", igraph::V(g_all)$local_id)

  nodes_df <- tibble::tibble(
    .name = igraph::V(g_all)$name,
    name = igraph::V(g_all)$name,
    family_id = igraph::V(g_all)$family_id,
    local_id = igraph::V(g_all)$local_id
  )

  families_df <- tibble::tibble(
    family_id = seq_along(family_sizes),
    size = family_sizes,
    type = purrr::map_chr(meta, "type")
  )

  family_objects <- families_df %>%
    dplyr::mutate(
      vertices = purrr::map(family_id, ~ igraph::V(g_all)$name[igraph::V(g_all)$family_id == .x]),
      g_sub = purrr::map(
        family_id,
        ~ igraph::induced_subgraph(
          g_all,
          vids = igraph::V(g_all)$name[igraph::V(g_all)$family_id == .x]
        )
      ),
      fam_df = purrr::map(family_id, ~ dplyr::filter(nodes_df, family_id == .x))
    )

  list(
    g_all = g_all,
    nodes_df = nodes_df,
    families_df = families_df,
    family_objects = family_objects
  )
}
