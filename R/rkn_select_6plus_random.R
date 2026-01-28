rkn_select_6plus_random <- function(fam_df,
                                    g_sub,
                                    seed = NULL,
                                    return_diagnostics = FALSE) {

  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  if (is.null(g_sub)) stop("g_sub non può essere NULL.")
  if (is.null(igraph::V(g_sub)$name)) stop("g_sub deve avere igraph::V(g_sub)$name valorizzato.")
  if (!all(fam_df$.name %in% igraph::V(g_sub)$name)) {
    stop("fam_df$.name contiene nodi non presenti in g_sub.")
  }

  n0 <- nrow(fam_df)
  if (n0 < 6) {
    warning("rkn_select_6plus_random chiamata con n < 6. Restituisco fallback: tutti discarded.")
    kept_final <- character(0)
    discarded_final <- fam_df$.name
    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))
    return(list(
      structure_name = "n_lt_6_fallback",
      kept_apriori = character(0),
      discarded_apriori = character(0),
      simple_candidates = character(0),
      complex_candidates = character(0),
      unknown = fam_df$.name,
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = FALSE
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  all_names <- fam_df$.name
  kept_apriori <- character(0)
  discarded_apriori <- character(0)

  discard_neighbors_of_kept <- function(g, kept_vec, discarded_vec) {
    if (length(kept_vec) == 0) return(discarded_vec)
    nb_idx <- igraph::neighborhood(g, order = 1, nodes = kept_vec)
    nb_names <- unique(unlist(lapply(nb_idx, function(vs) igraph::V(g)$name[vs])))
    nb_names <- setdiff(nb_names, kept_vec)
    unique(c(discarded_vec, nb_names))
  }

  # -------- 1) Riduzioni deterministiche (loop)
  changed <- TRUE
  while (changed) {
    changed <- FALSE

    # Propagazione: kept => discard neighbors
    new_discarded <- discard_neighbors_of_kept(g_sub, kept_apriori, discarded_apriori)
    if (length(setdiff(new_discarded, discarded_apriori)) > 0) {
      discarded_apriori <- new_discarded
      changed <- TRUE
    }

    undecided <- setdiff(all_names, union(kept_apriori, discarded_apriori))
    if (length(undecided) == 0) break

    g_u <- igraph::induced_subgraph(g_sub, vids = undecided)
    deg <- igraph::degree(g_u)

    # Isolati => kept
    iso <- names(deg)[deg == 0]
    if (length(iso) > 0) {
      kept_apriori <- unique(c(kept_apriori, iso))
      changed <- TRUE
      next
    }

    # Leaf forcing: leaf => keep(leaf), discard(neighbor)
    leaves <- names(deg)[deg == 1]
    if (length(leaves) > 0) {
      for (lf in leaves) {
        if (!(lf %in% igraph::V(g_u)$name)) next
        nb <- igraph::neighbors(g_u, v = lf, mode = "all")
        nb_name <- igraph::V(g_u)$name[nb][1]

        if (nb_name %in% kept_apriori) {
          discarded_apriori <- unique(c(discarded_apriori, lf))
        } else {
          kept_apriori <- unique(c(kept_apriori, lf))
          discarded_apriori <- unique(c(discarded_apriori, nb_name))
        }
        changed <- TRUE
      }
      next
    }

    # Star-hub puro: hub deg>=3 e tutti i vicini sono foglie (deg=1)
    if (igraph::vcount(g_u) >= 4) {
      deg_u <- igraph::degree(g_u)
      hubs <- names(deg_u)[deg_u >= 3]
      if (length(hubs) > 0) {
        for (h in hubs) {
          nbs <- igraph::neighbors(g_u, h)
          nb_names <- igraph::V(g_u)$name[nbs]
          if (length(nb_names) >= 3) {
            nb_deg <- deg_u[nb_names]
            if (all(nb_deg == 1)) {
              discarded_apriori <- unique(c(discarded_apriori, h))
              kept_apriori <- unique(c(kept_apriori, nb_names))
              changed <- TRUE
            }
          }
        }
        if (changed) next
      }
    }
  }

  # Propagazione finale
  discarded_apriori <- discard_neighbors_of_kept(g_sub, kept_apriori, discarded_apriori)
  undecided <- setdiff(all_names, union(kept_apriori, discarded_apriori))

  # Se non resta nulla
  if (length(undecided) == 0) {
    kept_final <- unique(kept_apriori)
    discarded_final <- setdiff(unique(discarded_apriori), kept_final)
    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))
    return(list(
      structure_name = "reduced_to_empty",
      kept_apriori = unique(kept_apriori),
      discarded_apriori = unique(discarded_apriori),
      simple_candidates = character(0),
      complex_candidates = character(0),
      unknown = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = FALSE
    ))
  }

  # -------- 2) Riconoscimento sul residuo
  g_r <- igraph::induced_subgraph(g_sub, vids = undecided)
  nr <- igraph::vcount(g_r)
  er <- igraph::ecount(g_r)
  degr <- igraph::degree(g_r)

  structure_name <- "unrecognized_6plus"
  used_largest_ivs <- FALSE

  simple_candidates <- character(0)
  complex_candidates <- character(0)
  unknown <- character(0)

  kept_afterrandomization <- character(0)
  discarded_afterrandomization <- character(0)

  pick_one_bipartite_side <- function(g) {
    b <- igraph::bipartite_mapping(g)
    if (!isTRUE(b$res)) return(character(0))
    A <- igraph::V(g)$name[b$type]
    B <- igraph::V(g)$name[!b$type]
    if (runif(1) < 0.5) A else B
  }

  # (i) star_n
  if (nr >= 6 && sum(degr == 1) == (nr - 1) && sum(degr == (nr - 1)) == 1) {
    structure_name <- paste0("star", nr)
    hub <- names(degr)[degr == (nr - 1)][1]
    leaves <- setdiff(igraph::V(g_r)$name, hub)
    kept_afterrandomization <- leaves
    discarded_afterrandomization <- hub

    # (ii) clique_n
  } else if (all(degr == (nr - 1))) {
    structure_name <- paste0("clique", nr)
    grp <- igraph::V(g_r)$name
    simple_candidates <- grp
    kept_afterrandomization <- sample(grp, 1)
    discarded_afterrandomization <- setdiff(grp, kept_afterrandomization)

    # (iii) path_n
  } else if (igraph::is_connected(g_r) &&
             er == (nr - 1) &&
             sum(degr == 1) == 2 &&
             all(degr %in% c(1, 2))) {

    structure_name <- paste0("path", nr)

    if (nr %% 2 == 1) {
      ends <- names(degr)[degr == 1]
      start <- ends[1]
      d <- igraph::distances(g_r, v = start)
      ord <- names(sort(d[1, ]))
      keep_idx <- seq(1, length(ord), by = 2)
      kept_afterrandomization <- ord[keep_idx]
      discarded_afterrandomization <- setdiff(ord, kept_afterrandomization)
    } else {
      complex_candidates <- igraph::V(g_r)$name
      chosen <- pick_one_bipartite_side(g_r)
      if (length(chosen) == 0) {
        unknown <- igraph::V(g_r)$name
      } else {
        kept_afterrandomization <- chosen
        discarded_afterrandomization <- setdiff(igraph::V(g_r)$name, kept_afterrandomization)
      }
    }

    # (iv) cycle_n
  } else if (igraph::is_connected(g_r) &&
             er == nr &&
             all(degr == 2)) {

    structure_name <- paste0("cycle", nr)

    if (nr %% 2 == 0) {
      complex_candidates <- igraph::V(g_r)$name
      chosen <- pick_one_bipartite_side(g_r)
      if (length(chosen) == 0) {
        unknown <- igraph::V(g_r)$name
      } else {
        kept_afterrandomization <- chosen
        discarded_afterrandomization <- setdiff(igraph::V(g_r)$name, kept_afterrandomization)
      }
    } else {
      unknown <- igraph::V(g_r)$name
    }

  } else {
    unknown <- igraph::V(g_r)$name
  }

  # -------- 3) Fallback MIS su unknown
  if (length(unknown) > 0) {
    used_largest_ivs <- TRUE
    g_ind <- igraph::induced_subgraph(g_sub, vids = unknown)
    sets <- igraph::largest_ivs(g_ind)
    set_names <- lapply(sets, function(s) igraph::V(g_ind)$name[s])
    chosen_names <- if (length(set_names) > 0) set_names[[sample(seq_along(set_names), 1)]] else character(0)

    kept_afterrandomization <- unique(c(kept_afterrandomization, chosen_names))
    discarded_afterrandomization <- unique(c(discarded_afterrandomization,
                                             setdiff(unknown, chosen_names)))
    complex_candidates <- unique(c(complex_candidates, unknown))
  }

  # -------- 4) Merge finali
  kept_final <- unique(c(kept_apriori, kept_afterrandomization))
  discarded_final <- unique(c(discarded_apriori, discarded_afterrandomization))
  discarded_final <- setdiff(discarded_final, kept_final)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  } else {
    return(list(
      structure_name = structure_name,
      kept_apriori = unique(kept_apriori),
      discarded_apriori = unique(discarded_apriori),
      simple_candidates = unique(simple_candidates),
      complex_candidates = unique(complex_candidates),
      unknown = unique(unknown),
      kept_afterrandomization = unique(kept_afterrandomization),
      discarded_afterrandomization = unique(discarded_afterrandomization),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = used_largest_ivs
    ))
  }
}

apply_rkn_select_6plus_random <- function(synth, seed = 999) {

  res <- synth$family_objects %>%
    dplyr::mutate(
      out = purrr::pmap(
        list(fam_df, g_sub, family_id),
        function(fam_df, g_sub, family_id) {
          rkn_select_6plus_random(
            fam_df = fam_df,
            g_sub = g_sub,
            seed = seed + as.integer(family_id),
            return_diagnostics = TRUE
          )
        }
      ),
      structure_name   = purrr::map_chr(out, "structure_name"),
      used_largest_ivs = purrr::map_lgl(out, "used_largest_ivs"),
      n_kept           = purrr::map_int(out, ~ length(.x$kept)),
      n_discarded      = purrr::map_int(out, ~ length(.x$discarded)),
      n_kept_apriori   = purrr::map_int(out, ~ length(.x$kept_apriori)),
      n_discarded_apriori = purrr::map_int(out, ~ length(.x$discarded_apriori)),
      n_simple         = purrr::map_int(out, ~ length(.x$simple_candidates)),
      n_complex        = purrr::map_int(out, ~ length(.x$complex_candidates)),
      n_unknown        = purrr::map_int(out, ~ length(.x$unknown))
    )

  summary_df <- res %>%
    dplyr::select(
      family_id, size, type, structure_name, used_largest_ivs,
      n_kept, n_discarded, n_kept_apriori, n_discarded_apriori,
      n_simple, n_complex, n_unknown
    )

  long_nodes <- res %>%
    dplyr::transmute(
      family_id,
      size,
      type,
      structure_name,
      used_largest_ivs,
      kept_apriori = purrr::map(out, "kept_apriori"),
      discarded_apriori = purrr::map(out, "discarded_apriori"),
      simple_candidates = purrr::map(out, "simple_candidates"),
      complex_candidates = purrr::map(out, "complex_candidates"),
      unknown = purrr::map(out, "unknown"),
      kept = purrr::map(out, "kept"),
      discarded = purrr::map(out, "discarded")
    ) %>%
    tidyr::pivot_longer(
      cols = c(kept_apriori, discarded_apriori, simple_candidates, complex_candidates, unknown, kept, discarded),
      names_to = "class",
      values_to = "nodes"
    ) %>%
    tidyr::unnest(nodes)

  list(
    results_by_family = res,
    summary_df = summary_df,
    long_nodes_df = long_nodes
  )
}
