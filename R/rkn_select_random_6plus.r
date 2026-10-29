#' Random selection for connected components with 6+ nodes
#'
#' @description
#' Applies the RedKinNet randomized node selection strategy to a connected
#' component with at least 6 nodes. The function combines:
#' \itemize{
#'   \item deterministic (a priori) decisions for recognizable structures; and
#'   \item randomized resolution for remaining ambiguous cases.
#' }
#'
#' The function expects the component to be provided both as a node-level
#' data frame and as an \pkg{igraph} subgraph; vertex names in the graph must
#' match node identifiers in the data frame.
#'
#' @param fam_df A data frame describing the nodes in the component. Must contain
#'   a character column \code{.name}. If \code{name} exists and \code{.name} does
#'   not, it is promoted to \code{.name}.
#' @param g_sub An \pkg{igraph} object representing the same component as
#'   \code{fam_df}. Must have \code{igraph::V(g_sub)$name}.
#' @param seed Optional integer RNG seed used for randomized decisions.
#' @param return_diagnostics Logical; if \code{TRUE}, returns additional
#'   per-node and per-family diagnostics.
#'
#' @return
#' If \code{return_diagnostics = FALSE}, a list with:
#' \itemize{
#'   \item \code{kept}: character vector of selected node names;
#'   \item \code{discarded}: character vector of discarded node names.
#' }
#'
#' If \code{return_diagnostics = TRUE}, the list also includes:
#' \itemize{
#'   \item \code{structure_name}: recognized structure label;
#'   \item \code{kept_apriori}, \code{discarded_apriori};
#'   \item \code{candidate_simple}, \code{candidate_complex};
#'   \item \code{kept_afterrandomization}, \code{discarded_afterrandomization};
#'   \item \code{mis_method};
#'   \item \code{node_diagnostics}: a data frame with one row per node;
#'   \item \code{family_diagnostics}: a one-row summary data frame.
#' }
#'
#' @details
#' This function is intended for internal pipeline use with controlled inputs.
#' It relies on \pkg{igraph} for graph operations and uses randomized tie-breaking
#' in specific structures. Set \code{seed} for reproducibility.
#'
#' @examples
#' \dontrun{
#' out <- rkn_select_random_6plus(fam_df, g_sub, seed = 1)
#' out_diag <- rkn_select_random_6plus(fam_df, g_sub, seed = 1, return_diagnostics = TRUE)
#' }
#'
#' @keywords internal
#' @noRd
rkn_select_random_6plus <- function(fam_df,
                                    g_sub,
                                    seed = NULL,
                                    return_diagnostics = FALSE) {

  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  if (is.null(g_sub)) stop("`g_sub` must not be NULL.")
  if (is.null(igraph::V(g_sub)$name)) stop("`g_sub` must have vertex names (igraph::V(g_sub)$name).")
  if (!all(fam_df$.name %in% igraph::V(g_sub)$name)) {
    stop("`fam_df$.name` contains nodes that are not in `g_sub`.")
  }

  n0 <- nrow(fam_df)
  if (n0 < 6) {
    warning("rkn_select_random_6plus() called on a component with fewer than 6 nodes; ",
            "all nodes are returned as discarded.")
    kept_final <- character(0)
    discarded_final <- fam_df$.name
    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    node_diag <- .rkn_df(
      id = fam_df$.name,
      pre_class = rep("discarded_apriori", length(fam_df$.name)),
      post_class = rep("discarded", length(fam_df$.name)),
      is_selected = rep(FALSE, length(fam_df$.name)),
      pre_reason = rep("n_lt_6_fallback", length(fam_df$.name))
    )
    fam_diag <- .rkn_df(
      n_nodes = n0,
      mis_method = NA_character_,
      structure_name = "n_lt_6_fallback"
    )

    return(list(
      structure_name = "n_lt_6_fallback",
      kept_apriori = character(0),
      discarded_apriori = fam_df$.name,
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      mis_method = NA_character_,
      node_diagnostics = node_diag,
      family_diagnostics = fam_diag
    ))
  }

  if (!is.null(seed)) {
    withr::local_preserve_seed()
    set.seed(seed)
  }

  all_names <- fam_df$.name
  kept_apriori <- character(0)
  discarded_apriori <- character(0)

  pre_reason_map <- stats::setNames(rep(NA_character_, length(all_names)), all_names)

  discard_neighbors_of_kept <- function(g, kept_vec, discarded_vec) {
    if (length(kept_vec) == 0) return(discarded_vec)
    nb_idx <- igraph::neighborhood(g, order = 1, nodes = kept_vec)
    g_names <- igraph::V(g)$name
    nb_names <- character(0)
    for (vs in nb_idx) nb_names <- c(nb_names, g_names[vs])
    nb_names <- setdiff(unique(nb_names), kept_vec)
    unique(c(discarded_vec, nb_names))
  }

  deg0 <- igraph::degree(g_sub, v = all_names)
  n_all <- length(all_names)

  if (n_all >= 6 &&
      sum(deg0 == 1) == (n_all - 1) &&
      sum(deg0 == (n_all - 1)) == 1) {

    structure_name <- paste0("star", n_all)
    hub <- names(deg0)[deg0 == (n_all - 1)][1]
    leaves <- setdiff(all_names, hub)

    kept_apriori <- leaves
    discarded_apriori <- hub

    for (x in leaves) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "star_deterministic_leaf"
    if (is.na(pre_reason_map[[hub]])) pre_reason_map[[hub]] <- "star_deterministic_hub"

    kept_final <- kept_apriori
    discarded_final <- discarded_apriori

    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    node_diag <- .rkn_build_node_diag_6plus(
      all_names = all_names,
      kept_apriori = kept_apriori,
      discarded_apriori = discarded_apriori,
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_final = kept_final,
      discarded_final = discarded_final,
      pre_reason_map = pre_reason_map
    )

    fam_diag <- .rkn_df(
      n_nodes = n_all,
      mis_method = NA_character_,
      structure_name = structure_name
    )

    return(list(
      structure_name = structure_name,
      kept_apriori = kept_apriori,
      discarded_apriori = discarded_apriori,
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      mis_method = NA_character_,
      node_diagnostics = node_diag,
      family_diagnostics = fam_diag
    ))
  }


  deg0 <- igraph::degree(g_sub, v = all_names)
  n_all <- length(all_names)

  if (n_all >= 6 &&
      igraph::is_connected(g_sub) &&
      igraph::ecount(g_sub) == (n_all - 1) &&
      sum(deg0 == 1) == 2 &&
      all(deg0 %in% c(1, 2))) {

    structure_name <- paste0("path", n_all)

    if (!is.null(seed)) set.seed(seed)

    ends <- names(deg0)[deg0 == 1]
    start <- ends[1]
    d <- igraph::distances(g_sub, v = start)
    ord <- names(sort(d[1, ]))  # ordine lungo il path

    if (n_all %% 2 == 1) {
      keep_idx <- seq(1, length(ord), by = 2)
      keep_set <- ord[keep_idx]
      disc_set <- setdiff(ord, keep_set)

      kept_apriori <- unique(c(kept_apriori, keep_set))
      discarded_apriori <- unique(c(discarded_apriori, disc_set))
      for (x in keep_set) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "path_odd_deterministic_keep"
      for (x in disc_set) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "path_odd_deterministic_discard"

      kept_final <- kept_apriori
      discarded_final <- discarded_apriori

      if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

      node_diag <- .rkn_build_node_diag_6plus(
        all_names = all_names,
        kept_apriori = kept_apriori,
        discarded_apriori = discarded_apriori,
        candidate_simple = character(0),
        candidate_complex = character(0),
        kept_final = kept_final,
        discarded_final = discarded_final,
        pre_reason_map = pre_reason_map
      )

      fam_diag <- .rkn_df(
        n_nodes = n_all,
        mis_method = NA_character_,
        structure_name = structure_name
      )

      return(list(
        structure_name = structure_name,
        kept_apriori = unique(kept_apriori),
        discarded_apriori = unique(discarded_apriori),
        candidate_simple = character(0),
        candidate_complex = character(0),
        kept_afterrandomization = character(0),
        discarded_afterrandomization = character(0),
        kept = kept_final,
        discarded = discarded_final,
        mis_method = NA_character_,
        node_diagnostics = node_diag,
        family_diagnostics = fam_diag
      ))

    } else {
      candidate_complex <- all_names
      for (x in candidate_complex) pre_reason_map[[x]] <- "path_even_bipartite"

      A <- ord[seq(1, length(ord), by = 2)]
      B <- ord[seq(2, length(ord), by = 2)]
      kept_afterrandomization <- if (stats::runif(1) < 0.5) A else B
      discarded_afterrandomization <- setdiff(all_names, kept_afterrandomization)

      kept_final <- kept_afterrandomization
      discarded_final <- discarded_afterrandomization

      if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

      node_diag <- .rkn_build_node_diag_6plus(
        all_names = all_names,
        kept_apriori = character(0),
        discarded_apriori = character(0),
        candidate_simple = character(0),
        candidate_complex = candidate_complex,
        kept_final = kept_final,
        discarded_final = discarded_final,
        pre_reason_map = pre_reason_map
      )

      fam_diag <- .rkn_df(
        n_nodes = n_all,
        mis_method = NA_character_,
        structure_name = structure_name
      )

      return(list(
        structure_name = structure_name,
        kept_apriori = character(0),
        discarded_apriori = character(0),
        candidate_simple = character(0),
        candidate_complex = candidate_complex,
        kept_afterrandomization = kept_final,
        discarded_afterrandomization = discarded_final,
        kept = kept_final,
        discarded = discarded_final,
        mis_method = NA_character_,
        node_diagnostics = node_diag,
        family_diagnostics = fam_diag
      ))
    }
  }

  deg0 <- igraph::degree(g_sub, v = all_names)
  n_all <- length(all_names)

  if (n_all >= 6 && all(deg0 == (n_all - 1))) {

    if (!is.null(seed)) set.seed(seed)

    structure_name <- paste0("clique", n_all)

    candidate_simple <- all_names
    candidate_complex <- character(0)

    for (x in candidate_simple) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "clique_choose1"

    kept_afterrandomization <- sample(candidate_simple, 1)
    discarded_afterrandomization <- setdiff(candidate_simple, kept_afterrandomization)

    kept_final <- kept_afterrandomization
    discarded_final <- discarded_afterrandomization

    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    node_diag <- .rkn_build_node_diag_6plus(
      all_names = all_names,
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidate_simple = candidate_simple,
      candidate_complex = character(0),
      kept_final = kept_final,
      discarded_final = discarded_final,
      pre_reason_map = pre_reason_map
    )

    fam_diag <- .rkn_df(
      n_nodes = n_all,
      mis_method = NA_character_,
      structure_name = structure_name
    )

    return(list(
      structure_name = structure_name,
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidate_simple = candidate_simple,
      candidate_complex = character(0),
      kept_afterrandomization = kept_afterrandomization,
      discarded_afterrandomization = discarded_afterrandomization,
      kept = kept_final,
      discarded = discarded_final,
      mis_method = NA_character_,
      node_diagnostics = node_diag,
      family_diagnostics = fam_diag
    ))
  }

  if (n_all >= 6 &&
      igraph::is_connected(g_sub) &&
      igraph::ecount(g_sub) == n_all &&
      all(deg0 == 2)) {

    if (!is.null(seed)) set.seed(seed)

    structure_name <- paste0("cycle", n_all)

    candidate_simple <- character(0)
    candidate_complex <- all_names

    if (n_all %% 2 == 0) {
      for (x in candidate_complex) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "cycle_even_bipartite"

      # Explicit cyclic ordering by traversal: ordering by graph distance is
      # insufficient for cycles, because several vertices share the same distance.
      start <- all_names[1]
      nbrs <- igraph::neighbors(g_sub, start, mode = "all")
      if (length(nbrs) != 2) stop("Cycle traversal invariant violated: start node degree is not 2.")

      prev <- start
      curr <- igraph::V(g_sub)$name[nbrs][1]  # choose an arbitrary direction

      ord <- character(n_all)
      ord[1] <- start

      for (k in 2:n_all) {
        ord[k] <- curr
        nb <- igraph::neighbors(g_sub, curr, mode = "all")
        nb_names <- igraph::V(g_sub)$name[nb]

        # Each vertex of a simple cycle has two neighbours: advance to the one not just visited.
        nxt <- setdiff(nb_names, prev)
        if (length(nxt) != 1) stop("Cycle traversal failed: non-unique next vertex.")
        prev <- curr
        curr <- nxt
      }

      # Alternating along the cyclic ordering guarantees independence.
      A <- ord[seq(1, n_all, by = 2)]
      B <- ord[seq(2, n_all, by = 2)]
      kept_afterrandomization <- if (stats::runif(1) < 0.5) A else B
      discarded_afterrandomization <- setdiff(all_names, kept_afterrandomization)

      mis_method <- NA_character_

    } else {
      for (x in candidate_complex) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "exact_fallback"

      drawn <- .rkn_draw_maximum_set(igraph::induced_subgraph(g_sub, vids = all_names))
      mis_method <- drawn$method

      kept_afterrandomization <- drawn$kept
      discarded_afterrandomization <- setdiff(all_names, kept_afterrandomization)
    }

    kept_final <- kept_afterrandomization
    discarded_final <- discarded_afterrandomization

    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    node_diag <- .rkn_build_node_diag_6plus(
      all_names = all_names,
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidate_simple = character(0),
      candidate_complex = candidate_complex,
      kept_final = kept_final,
      discarded_final = discarded_final,
      pre_reason_map = pre_reason_map
    )

    fam_diag <- .rkn_df(
      n_nodes = n_all,
      mis_method = mis_method,
      structure_name = structure_name
    )

    return(list(
      structure_name = structure_name,
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidate_simple = character(0),
      candidate_complex = candidate_complex,
      kept_afterrandomization = kept_afterrandomization,
      discarded_afterrandomization = discarded_afterrandomization,
      kept = kept_final,
      discarded = discarded_final,
      mis_method = mis_method,
      node_diagnostics = node_diag,
      family_diagnostics = fam_diag
    ))
  }

  changed <- TRUE
  while (changed) {
    changed <- FALSE

    new_discarded <- discard_neighbors_of_kept(g_sub, kept_apriori, discarded_apriori)
    if (length(setdiff(new_discarded, discarded_apriori)) > 0) {
      discarded_apriori <- new_discarded
      for (x in setdiff(new_discarded, names(pre_reason_map)[!is.na(pre_reason_map)])) {
        if (x %in% names(pre_reason_map) && is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "neighbor_of_kept"
      }
      changed <- TRUE
    }

    undecided <- setdiff(all_names, union(kept_apriori, discarded_apriori))
    if (length(undecided) == 0) break

    g_u <- igraph::induced_subgraph(g_sub, vids = undecided)
    deg <- igraph::degree(g_u)

    iso <- names(deg)[deg == 0]
    if (length(iso) > 0) {
      kept_apriori <- unique(c(kept_apriori, iso))
      for (x in iso) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "isolated"
      changed <- TRUE
      next
    }

    leaves <- names(deg)[deg == 1]
    if (length(leaves) > 0) {
      for (lf in leaves) {
        if (!(lf %in% igraph::V(g_u)$name)) next
        nb <- igraph::neighbors(g_u, v = lf, mode = "all")
        nb_name <- igraph::V(g_u)$name[nb][1]

        if (nb_name %in% kept_apriori) {
          discarded_apriori <- unique(c(discarded_apriori, lf))
          if (is.na(pre_reason_map[[lf]])) pre_reason_map[[lf]] <- "leaf_neighbor_is_kept"
        } else {
          kept_apriori <- unique(c(kept_apriori, lf))
          discarded_apriori <- unique(c(discarded_apriori, nb_name))
          if (is.na(pre_reason_map[[lf]])) pre_reason_map[[lf]] <- "leaf_forcing_keep_leaf"
          if (is.na(pre_reason_map[[nb_name]])) pre_reason_map[[nb_name]] <- "leaf_forcing_discard_neighbor"
        }
        changed <- TRUE
      }
      next
    }

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
              if (is.na(pre_reason_map[[h]])) pre_reason_map[[h]] <- "pure_star_hub"
              for (x in nb_names) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "pure_star_leaf"
              changed <- TRUE
            }
          }
        }
        if (changed) next
      }
    }
  }

  discarded_apriori <- discard_neighbors_of_kept(g_sub, kept_apriori, discarded_apriori)
  undecided <- setdiff(all_names, union(kept_apriori, discarded_apriori))

  if (length(undecided) == 0) {
    kept_final <- unique(kept_apriori)
    discarded_final <- setdiff(unique(discarded_apriori), kept_final)
    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    node_diag <- .rkn_build_node_diag_6plus(
      all_names = all_names,
      kept_apriori = kept_apriori,
      discarded_apriori = discarded_apriori,
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_final = kept_final,
      discarded_final = discarded_final,
      pre_reason_map = pre_reason_map
    )
    fam_diag <- .rkn_df(
      n_nodes = n0,
      mis_method = NA_character_,
      structure_name = "reduced_to_empty"
    )

    return(list(
      structure_name = "reduced_to_empty",
      kept_apriori = unique(kept_apriori),
      discarded_apriori = unique(discarded_apriori),
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      mis_method = NA_character_,
      node_diagnostics = node_diag,
      family_diagnostics = fam_diag
    ))
  }

  g_r <- igraph::induced_subgraph(g_sub, vids = undecided)
  nr <- igraph::vcount(g_r)
  er <- igraph::ecount(g_r)
  degr <- igraph::degree(g_r)

  structure_name <- "unrecognized_6plus"
  mis_method <- NA_character_

  candidate_simple <- character(0)
  candidate_complex <- character(0)

  fallback_nodes <- character(0)

  kept_afterrandomization <- character(0)
  discarded_afterrandomization <- character(0)

  pick_one_bipartite_side <- function(g) {
    b <- igraph::bipartite_mapping(g)
    if (!isTRUE(b$res)) return(character(0))
    A <- igraph::V(g)$name[b$type]
    B <- igraph::V(g)$name[!b$type]
    if (stats::runif(1) < 0.5) A else B
  }

  if (nr >= 6 && sum(degr == 1) == (nr - 1) && sum(degr == (nr - 1)) == 1) {
    structure_name <- paste0("star", nr)
    hub <- names(degr)[degr == (nr - 1)][1]
    leaves <- setdiff(igraph::V(g_r)$name, hub)

    kept_apriori <- unique(c(kept_apriori, leaves))
    discarded_apriori <- unique(c(discarded_apriori, hub))
    for (x in leaves) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "star_deterministic_leaf"
    if (is.na(pre_reason_map[[hub]])) pre_reason_map[[hub]] <- "star_deterministic_hub"

  } else if (all(degr == (nr - 1))) {
    structure_name <- paste0("clique", nr)
    grp <- igraph::V(g_r)$name
    candidate_simple <- grp
    for (x in grp) pre_reason_map[[x]] <- "clique_choose1"

    kept_afterrandomization <- sample(grp, 1)
    discarded_afterrandomization <- setdiff(grp, kept_afterrandomization)

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
      keep_set <- ord[keep_idx]
      disc_set <- setdiff(ord, keep_set)

      kept_apriori <- unique(c(kept_apriori, keep_set))
      discarded_apriori <- unique(c(discarded_apriori, disc_set))
      for (x in keep_set) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "path_odd_deterministic_keep"
      for (x in disc_set) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "path_odd_deterministic_discard"
    } else {
      candidate_complex <- igraph::V(g_r)$name
      for (x in candidate_complex) pre_reason_map[[x]] <- "path_even_bipartite"

      chosen <- pick_one_bipartite_side(g_r)
      if (length(chosen) == 0) {
        fallback_nodes <- igraph::V(g_r)$name
      } else {
        kept_afterrandomization <- chosen
        discarded_afterrandomization <- setdiff(igraph::V(g_r)$name, kept_afterrandomization)
      }
    }

  } else if (igraph::is_connected(g_r) &&
             er == nr &&
             all(degr == 2)) {

    structure_name <- paste0("cycle", nr)

    if (nr %% 2 == 0) {
      candidate_complex <- igraph::V(g_r)$name
      for (x in candidate_complex) pre_reason_map[[x]] <- "cycle_even_bipartite"

      chosen <- pick_one_bipartite_side(g_r)
      if (length(chosen) == 0) {
        fallback_nodes <- igraph::V(g_r)$name
      } else {
        kept_afterrandomization <- chosen
        discarded_afterrandomization <- setdiff(igraph::V(g_r)$name, kept_afterrandomization)
      }
    } else {
      fallback_nodes <- igraph::V(g_r)$name
    }

  } else {
    fallback_nodes <- igraph::V(g_r)$name
  }

  if (length(fallback_nodes) > 0) {
    candidate_complex <- unique(c(candidate_complex, fallback_nodes))
    for (x in fallback_nodes) pre_reason_map[[x]] <- "exact_fallback"

    drawn <- .rkn_draw_maximum_set(igraph::induced_subgraph(g_sub, vids = fallback_nodes))
    mis_method <- drawn$method
    chosen_names <- drawn$kept

    kept_afterrandomization <- unique(c(kept_afterrandomization, chosen_names))
    discarded_afterrandomization <- unique(c(discarded_afterrandomization,
                                             setdiff(fallback_nodes, chosen_names)))
  }

  kept_final <- unique(c(kept_apriori, kept_afterrandomization))
  discarded_final <- unique(c(discarded_apriori, discarded_afterrandomization))
  discarded_final <- setdiff(discarded_final, kept_final)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  }

  node_diag <- .rkn_build_node_diag_6plus(
    all_names = all_names,
    kept_apriori = kept_apriori,
    discarded_apriori = discarded_apriori,
    candidate_simple = candidate_simple,
    candidate_complex = candidate_complex,
    kept_final = kept_final,
    discarded_final = discarded_final,
    pre_reason_map = pre_reason_map
  )

  fam_diag <- .rkn_df(
    n_nodes = n0,
    mis_method = mis_method,
    structure_name = structure_name
  )

  return(list(
    structure_name = structure_name,
    kept_apriori = unique(kept_apriori),
    discarded_apriori = unique(discarded_apriori),
    candidate_simple = unique(candidate_simple),
    candidate_complex = unique(candidate_complex),
    kept_afterrandomization = unique(kept_afterrandomization),
    discarded_afterrandomization = unique(discarded_afterrandomization),
    kept = kept_final,
    discarded = discarded_final,
    mis_method = mis_method,
    node_diagnostics = node_diag,
    family_diagnostics = fam_diag
  ))
}

#' Internal helper: build per-node diagnostics for 6+ selection
#'
#' @description Constructs a per-node diagnostics data frame from intermediate and final classification vectors.
#'
#' @keywords internal
#' @noRd
.rkn_build_node_diag_6plus <- function(all_names,
                                       kept_apriori,
                                       discarded_apriori,
                                       candidate_simple,
                                       candidate_complex,
                                       kept_final,
                                       discarded_final,
                                       pre_reason_map = NULL) {
  ids <- all_names

  pre_class <- stats::setNames(rep("candidate_complex", length(ids)), ids)
  pre_class[candidate_simple] <- "candidate_simple"
  pre_class[candidate_complex] <- "candidate_complex"
  pre_class[discarded_apriori] <- "discarded_apriori"
  pre_class[kept_apriori] <- "kept_apriori"

  post_class <- stats::setNames(rep(NA_character_, length(ids)), ids)
  post_class[kept_final] <- "kept"
  post_class[discarded_final] <- "discarded"

  if (is.null(pre_reason_map)) {
    pre_reason <- rep(NA_character_, length(ids))
  } else {
    if (is.null(names(pre_reason_map))) names(pre_reason_map) <- ids
    pre_reason <- unname(pre_reason_map[ids])
  }

  .rkn_df(
    id = ids,
    pre_class = unname(pre_class[ids]),
    post_class = unname(post_class[ids]),
    is_selected = unname(post_class[ids] == "kept"),
    pre_reason = pre_reason
  )
}
