#' Internal helper: recognize a canonical structure and build its candidates (6+)
#'
#' @description
#' Classifies a component (or post-peeling residual) as one of the structures
#' \code{rkn_select_random_6plus()} recognizes -- star, clique, path, (bipartite or odd)
#' cycle -- or \code{"irregular"} if none apply. For fully deterministic structures
#' (star, odd path), returns the resolved \code{kept}/\code{discarded} directly. For
#' ambiguous structures (clique, even path, cycle), returns \code{candidate_sets}: for a
#' clique, one candidate per single node; for a bipartite path/cycle, its two sides
#' (found via \code{igraph::bipartite_mapping()}, which works uniformly for both since
#' both are bipartite); for an odd cycle, its maximum independent sets, listed by the
#' exact engine.
#'
#' @keywords internal
#' @noRd
.rkn_recognize_candidates_6plus <- function(g) {
  all_names <- igraph::V(g)$name
  n <- length(all_names)
  deg <- igraph::degree(g, v = all_names)

  if (n >= 2 && sum(deg == 1) == (n - 1) && sum(deg == (n - 1)) == 1) {
    hub <- names(deg)[deg == (n - 1)][1]
    leaves <- setdiff(all_names, hub)
    return(list(
      kind = "star", label = "star", n = n, candidate_sets = NULL,
      kept = leaves, discarded = hub,
      pre_reason = stats::setNames(c(rep("star_deterministic_leaf", length(leaves)), "star_deterministic_hub"),
                            c(leaves, hub))
    ))
  }

  if (n >= 2 && all(deg == (n - 1))) {
    return(list(
      kind = "clique", label = "clique", n = n, candidate_sets = as.list(all_names),
      kept = character(0), discarded = character(0),
      pre_reason = stats::setNames(rep("clique_choose1", n), all_names)
    ))
  }

  if (igraph::is_connected(g) && igraph::ecount(g) == (n - 1) &&
      sum(deg == 1) == 2 && all(deg %in% c(1, 2))) {

    if (n %% 2 == 1) {
      ends <- names(deg)[deg == 1]
      d <- igraph::distances(g, v = ends[1])
      ord <- names(sort(d[1, ]))
      keep_set <- ord[seq(1, n, by = 2)]
      disc_set <- setdiff(ord, keep_set)
      return(list(
        kind = "path_odd", label = "path", n = n, candidate_sets = NULL,
        kept = keep_set, discarded = disc_set,
        pre_reason = stats::setNames(c(rep("path_odd_deterministic_keep", length(keep_set)),
                                rep("path_odd_deterministic_discard", length(disc_set))),
                              c(keep_set, disc_set))
      ))
    }

    bp <- igraph::bipartite_mapping(g)
    A <- all_names[bp$type]; B <- all_names[!bp$type]
    return(list(
      kind = "path_even", label = "path", n = n, candidate_sets = list(A, B),
      kept = character(0), discarded = character(0),
      pre_reason = stats::setNames(rep("path_even_bipartite", n), all_names)
    ))
  }

  # A cycle is bipartite iff n is even; checked directly rather than via
  # igraph::bipartite_mapping()'s own res flag, which errors (rather than
  # returning res = FALSE) on a non-bipartite graph in some igraph versions.
  if (igraph::is_connected(g) && igraph::ecount(g) == n && all(deg == 2)) {
    if (n %% 2 == 0) {
      bp <- igraph::bipartite_mapping(g)
      A <- all_names[bp$type]; B <- all_names[!bp$type]
      return(list(
        kind = "cycle_even", label = "cycle", n = n, candidate_sets = list(A, B),
        kept = character(0), discarded = character(0),
        pre_reason = stats::setNames(rep("cycle_even_bipartite", n), all_names)
      ))
    }
    sets <- .rkn_mis_enumerate(.rkn_mis_exact(.rkn_adjacency(g)), limit = Inf)
    candidate_sets <- vector("list", ncol(sets))
    for (j in seq_len(ncol(sets))) candidate_sets[[j]] <- all_names[sets[, j]]
    return(list(
      kind = "cycle_odd", label = "cycle", n = n, candidate_sets = candidate_sets,
      kept = character(0), discarded = character(0),
      pre_reason = stats::setNames(rep("exact_cycle_odd", n), all_names)
    ))
  }

  list(kind = "irregular", label = NA_character_, n = n, candidate_sets = NULL,
      kept = character(0), discarded = character(0), pre_reason = character(0))
}

#' Criteria-based selection for connected components with 6+ nodes
#'
#' @description
#' Applies the same structural recognition and deterministic peeling logic as
#' \code{rkn_select_random_6plus()} to a connected component with at least 6 nodes,
#' but resolves every ambiguous case (an even path or cycle's two bipartite sides, a
#' clique's single choice, an odd cycle's maximum independent sets) by criterion score
#' instead of a uniform random draw, via the same \code{.rkn_select_among_candidates()}
#' engine used by \code{rkn_select_criteria_upto5()}. Ties are broken with a seeded
#' random draw restricted to the tied candidates. \code{structure_name} uses the same
#' coarse labels as \code{rkn_select_random_6plus()} (\code{"pathN"}, \code{"cycleN"}, not
#' distinguishing odd/even in the name).
#'
#' @param fam_df A data frame describing the nodes in the component. Must contain a
#'   character column \code{.name}. If \code{name} exists and \code{.name} does not, it
#'   is promoted to \code{.name}.
#' @param g_sub An \pkg{igraph} object representing the same component as \code{fam_df}.
#'   Must have \code{igraph::V(g_sub)$name}.
#' @param criterion_primary An \code{rkn_criterion} object (see
#'   \code{\link{rkn_criterion_target_mean}} and siblings) used to score candidate
#'   resolutions.
#' @param criterion_secondary Optional \code{rkn_criterion} object used to break ties in
#'   \code{criterion_primary}.
#' @param seed Optional integer RNG seed, used only to break exact ties in the criterion
#'   score(s).
#' @param return_diagnostics Logical; if \code{TRUE}, returns additional per-node and
#'   per-family diagnostics.
#'
#' @return Same structure as \code{rkn_select_criteria_upto5()}.
#'
#' @seealso \code{rkn_select_random_6plus()}, \code{rkn_select_criteria_upto5()}
#'
#' @keywords internal
#' @noRd
rkn_select_criteria_6plus <- function(fam_df,
                                      g_sub,
                                      criterion_primary,
                                      criterion_secondary = NULL,
                                      seed = NULL,
                                      return_diagnostics = FALSE) {

  if (missing(criterion_primary) || !inherits(criterion_primary, "rkn_criterion")) {
    stop("`criterion_primary` must be an rkn_criterion object (see rkn_criterion_*()).")
  }
  if (!is.null(criterion_secondary) && !inherits(criterion_secondary, "rkn_criterion")) {
    stop("`criterion_secondary` must be an rkn_criterion object or NULL.")
  }

  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  if (is.null(g_sub)) stop("`g_sub` must not be NULL.")
  if (is.null(igraph::V(g_sub)$name)) stop("`g_sub` must have `igraph::V(g_sub)$name` set.")
  if (!all(fam_df$.name %in% igraph::V(g_sub)$name)) {
    stop("`fam_df$.name` contains nodes not present in `g_sub`.")
  }

  all_names <- fam_df$.name
  n0 <- length(all_names)

  if (n0 < 6) {
    warning("rkn_select_criteria_6plus() called with n < 6. Returning fallback: all discarded.")
    kept_final <- character(0)
    discarded_final <- all_names
    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    node_diag <- .rkn_build_node_diag_6plus(
      all_names = all_names, kept_apriori = character(0), discarded_apriori = all_names,
      candidate_simple = character(0), candidate_complex = character(0),
      kept_final = kept_final, discarded_final = discarded_final,
      pre_reason_map = stats::setNames(rep("n_lt_6_fallback", n0), all_names)
    )
    fam_diag <- .rkn_df(n_nodes = n0, mis_method = NA_character_, structure_name = "n_lt_6_fallback")
    return(list(
      structure_name = "n_lt_6_fallback", kept_apriori = character(0), discarded_apriori = all_names,
      candidate_simple = character(0), candidate_complex = character(0),
      kept_aftercriteria = character(0), discarded_aftercriteria = character(0),
      kept = kept_final, discarded = discarded_final, mis_method = NA_character_,
      node_diagnostics = node_diag, family_diagnostics = fam_diag
    ))
  }

  # Resolve a graph (whole component, or post-peeling residual) against a recognized
  # structure. Returns NULL if the graph is irregular (caller handles the generic
  # exact-engine fallback, since only the residual case needs it). `fixed_kept` holds the
  # nodes already kept by the peeling loop, so that candidates are scored on the
  # family-level kept set, as in rkn_select_criteria_upto5().
  resolve_component <- function(g, seed_offset, fixed_kept = character(0)) {
    shape <- .rkn_recognize_candidates_6plus(g)
    if (shape$kind == "irregular") return(NULL)

    if (is.null(shape$candidate_sets)) {
      return(list(
        structure_name = paste0(shape$label, shape$n), mis_method = NA_character_,
        kept_apriori = shape$kept, discarded_apriori = shape$discarded,
        kept_aftercriteria = character(0), discarded_aftercriteria = character(0),
        candidate_simple = character(0), candidate_complex = character(0),
        pre_reason = shape$pre_reason
      ))
    }

    picked <- .rkn_select_among_candidates(
      candidate_sets = shape$candidate_sets, fam_df = fam_df, fixed_kept_ids = fixed_kept,
      criterion_primary = criterion_primary, criterion_secondary = criterion_secondary,
      seed = if (!is.null(seed)) seed + seed_offset else NULL
    )
    kept <- picked$chosen_set
    all_in_shape <- unique(unlist(shape$candidate_sets))
    discarded <- setdiff(all_in_shape, kept)
    is_complex <- shape$kind != "clique"

    list(
      structure_name = paste0(shape$label, shape$n),
      mis_method = if (identical(shape$kind, "cycle_odd")) "exact" else NA_character_,
      kept_apriori = character(0), discarded_apriori = character(0),
      kept_aftercriteria = kept, discarded_aftercriteria = discarded,
      candidate_simple = if (!is_complex) all_in_shape else character(0),
      candidate_complex = if (is_complex) all_in_shape else character(0),
      pre_reason = shape$pre_reason
    )
  }

  finalize <- function(res, pre_reason_extra = NULL) {
    kept_final <- unique(c(res$kept_apriori, res$kept_aftercriteria))
    discarded_final <- unique(c(res$discarded_apriori, res$discarded_aftercriteria))
    discarded_final <- setdiff(discarded_final, kept_final)

    if (!return_diagnostics) return(list(kept = kept_final, discarded = discarded_final))

    # Individuals in every or in no maximum set of the residual searched were decided
    # without a choice, whatever the criteria.
    if (!is.null(res$inclusion)) {
      prob <- res$inclusion
      res$kept_apriori <- unique(c(res$kept_apriori, names(prob)[prob == 1]))
      res$discarded_apriori <- unique(c(res$discarded_apriori, names(prob)[prob == 0]))
      res$candidate_complex <- setdiff(res$candidate_complex, names(prob)[prob %in% c(0, 1)])
    }

    pre_reason_map <- stats::setNames(rep(NA_character_, n0), all_names)
    if (!is.null(pre_reason_extra)) pre_reason_map[names(pre_reason_extra)] <- pre_reason_extra
    pre_reason_map[names(res$pre_reason)] <- res$pre_reason

    node_diag <- .rkn_build_node_diag_6plus(
      all_names = all_names, kept_apriori = res$kept_apriori, discarded_apriori = res$discarded_apriori,
      candidate_simple = res$candidate_simple, candidate_complex = res$candidate_complex,
      kept_final = kept_final, discarded_final = discarded_final, pre_reason_map = pre_reason_map
    )
    fam_diag <- .rkn_df(n_nodes = n0, mis_method = res$mis_method, structure_name = res$structure_name)

    list(
      structure_name = res$structure_name, kept_apriori = unique(res$kept_apriori),
      discarded_apriori = unique(res$discarded_apriori),
      candidate_simple = unique(res$candidate_simple), candidate_complex = unique(res$candidate_complex),
      kept_aftercriteria = unique(res$kept_aftercriteria), discarded_aftercriteria = unique(res$discarded_aftercriteria),
      kept = kept_final, discarded = discarded_final, mis_method = res$mis_method,
      node_diagnostics = node_diag, family_diagnostics = fam_diag
    )
  }

  g_whole <- igraph::induced_subgraph(g_sub, vids = all_names)
  res_whole <- resolve_component(g_whole, seed_offset = 0L)
  if (!is.null(res_whole)) return(finalize(res_whole))

  kept_apriori <- character(0)
  discarded_apriori <- character(0)
  pre_reason_map <- stats::setNames(rep(NA_character_, n0), all_names)

  exclude_neighbors <- function(kept_vec, discarded_vec) {
    if (length(kept_vec) == 0) return(discarded_vec)
    nb_idx <- igraph::neighborhood(g_sub, order = 1, nodes = kept_vec)
    sub_names <- igraph::V(g_sub)$name
    nb_names <- character(0)
    for (vs in nb_idx) nb_names <- c(nb_names, sub_names[vs])
    nb_names <- unique(setdiff(nb_names, kept_vec))
    unique(c(discarded_vec, nb_names))
  }

  changed <- TRUE
  while (changed) {
    changed <- FALSE

    new_discarded <- exclude_neighbors(kept_apriori, discarded_apriori)
    if (length(setdiff(new_discarded, discarded_apriori)) > 0) {
      for (x in setdiff(new_discarded, discarded_apriori)) {
        if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "neighbor_of_kept"
      }
      discarded_apriori <- new_discarded
      changed <- TRUE
    }

    undecided <- setdiff(all_names, union(kept_apriori, discarded_apriori))
    if (length(undecided) == 0) break

    g_u <- igraph::induced_subgraph(g_sub, vids = undecided)
    deg_u <- igraph::degree(g_u)

    iso <- names(deg_u)[deg_u == 0]
    if (length(iso) > 0) {
      kept_apriori <- unique(c(kept_apriori, iso))
      for (x in iso) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "isolated"
      changed <- TRUE
      next
    }

    leaves <- names(deg_u)[deg_u == 1]
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
      hubs <- names(deg_u)[deg_u >= 3]
      for (h in hubs) {
        nbs <- igraph::neighbors(g_u, h)
        nb_names <- igraph::V(g_u)$name[nbs]
        if (length(nb_names) >= 3 && all(deg_u[nb_names] == 1)) {
          discarded_apriori <- unique(c(discarded_apriori, h))
          kept_apriori <- unique(c(kept_apriori, nb_names))
          if (is.na(pre_reason_map[[h]])) pre_reason_map[[h]] <- "pure_star_hub"
          for (x in nb_names) if (is.na(pre_reason_map[[x]])) pre_reason_map[[x]] <- "pure_star_leaf"
          changed <- TRUE
        }
      }
    }
  }

  discarded_apriori <- exclude_neighbors(kept_apriori, discarded_apriori)
  undecided <- setdiff(all_names, union(kept_apriori, discarded_apriori))

  if (length(undecided) == 0) {
    return(finalize(list(
      structure_name = "reduced_to_empty", mis_method = NA_character_,
      kept_apriori = kept_apriori, discarded_apriori = discarded_apriori,
      kept_aftercriteria = character(0), discarded_aftercriteria = character(0),
      candidate_simple = character(0), candidate_complex = character(0), pre_reason = character(0)
    ), pre_reason_extra = pre_reason_map))
  }

  g_r <- igraph::induced_subgraph(g_sub, vids = undecided)
  res_r <- resolve_component(g_r, seed_offset = 1L, fixed_kept = kept_apriori)

  if (is.null(res_r)) {
    picked <- .rkn_criteria_maximum_set(
      g_r, fam_df = fam_df, fixed_kept = kept_apriori,
      criterion_primary = criterion_primary, criterion_secondary = criterion_secondary,
      seed = if (!is.null(seed)) seed + 2L else NULL
    )
    res_r <- list(
      structure_name = "unrecognized_6plus", mis_method = picked$method,
      kept_apriori = character(0), discarded_apriori = character(0),
      kept_aftercriteria = picked$kept,
      discarded_aftercriteria = setdiff(igraph::V(g_r)$name, picked$kept),
      candidate_simple = character(0), candidate_complex = igraph::V(g_r)$name,
      inclusion = if (return_diagnostics) .rkn_search_inclusion(g_r, picked$searched)$inclusion,
      pre_reason = stats::setNames(rep("exact_fallback", igraph::vcount(g_r)), igraph::V(g_r)$name)
    )
  }

  combined <- list(
    structure_name = res_r$structure_name, mis_method = res_r$mis_method,
    kept_apriori = unique(c(kept_apriori, res_r$kept_apriori)),
    discarded_apriori = unique(c(discarded_apriori, res_r$discarded_apriori)),
    kept_aftercriteria = res_r$kept_aftercriteria, discarded_aftercriteria = res_r$discarded_aftercriteria,
    candidate_simple = res_r$candidate_simple, candidate_complex = res_r$candidate_complex,
    pre_reason = res_r$pre_reason, inclusion = res_r$inclusion
  )
  finalize(combined, pre_reason_extra = pre_reason_map)
}
