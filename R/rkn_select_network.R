.rkn_reserved_attributes <- c("degree", "max_path_length", ".name")

#' Internal helper: error message for vertex attributes with reserved names
#'
#' The structural rules read the degree and the eccentricity of each individual from
#' columns named \code{degree} and \code{max_path_length} of the family table, next to
#' the vertex attributes; an attribute with one of these names would be mistaken for
#' them.
#'
#' @keywords internal
#' @noRd
.rkn_reserved_message <- function(reserved) {
  paste0("Vertex attributes cannot be named ", paste0("\"", .rkn_reserved_attributes, "\"", collapse = ", "),
         ", which are reserved for the structural measures computed during the selection; rename ",
         paste0("\"", reserved, "\"", collapse = ", "), ".")
}

#' Internal helper: network-level selection shared by the random and criteria selectors
#'
#' Partitions \code{g} into families (connected components), resolves each family with
#' \code{select_upto5} (families of at most 5 individuals) or \code{select_6plus}
#' (larger families), and returns the selection as an \code{rkn_selection} object, with
#' the diagnostics tables when requested. Both per-family functions are called with the
#' arguments \code{fam_df}, \code{g_sub}, \code{seed} and \code{return_diagnostics};
#' \code{method}, \code{criteria} and \code{call} are stored in the object. The
#' arguments are validated by the exported selectors.
#'
#' @keywords internal
#' @noRd
.rkn_select_network <- function(g, select_upto5, select_6plus, seed, diagnostics,
                                method, criteria, call, uniform = FALSE) {
  n_nodes <- igraph::vcount(g)

  # Families are extracted by integer vertex id: extraction by name scans the
  # whole graph for every family and makes the run time quadratic.
  fam_idx <- split(seq_len(n_nodes), igraph::components(g)$membership)
  n_fam <- length(fam_idx)
  fam_sizes <- lengths(fam_idx)

  # Extract a per-family data.frame of vertex attributes (robust to missing attributes).
  # Node attributes and the structural measures required by rules_db are computed once
  # on the whole graph: repeating these igraph calls for every family dominates the run
  # time on large networks. Degree and eccentricity (max_path_length) are within-family
  # quantities, since igraph ignores vertices in other components.
  node_cols <- c(list(.name = igraph::V(g)$name), igraph::vertex_attr(g))
  if (is.null(node_cols[["name"]])) node_cols[["name"]] <- node_cols[[".name"]]
  node_cols[["degree"]] <- as.integer(igraph::degree(g, mode = "all"))
  node_cols[["max_path_length"]] <- as.integer(igraph::eccentricity(g, mode = "all"))

  .fam_df_from_ids <- function(ids) {
    cols <- node_cols
    for (j in seq_along(cols)) cols[[j]] <- cols[[j]][ids]
    do.call(.rkn_df, cols)
  }

  kept_list <- vector("list", n_fam)
  discarded_list <- vector("list", n_fam)

  mis_method_vec <- rep(NA_character_, n_fam)
  uniform_vec <- rep(FALSE, n_fam)
  n_max_sets_vec <- rep(NA_real_, n_fam)
  node_diag_list <- if (diagnostics) vector("list", n_fam) else NULL
  if (diagnostics) {
    fs_selector <- character(n_fam)
    fs_structure <- rep(NA_character_, n_fam)
    fs_n_kept <- integer(n_fam)
    fs_n_discarded <- integer(n_fam)
    fs_n_kept_apriori <- rep(NA_integer_, n_fam)
    fs_n_discarded_apriori <- rep(NA_integer_, n_fam)
    fs_n_candidate_simple <- rep(NA_integer_, n_fam)
    fs_n_candidate_complex <- rep(NA_integer_, n_fam)
  }

  for (i in seq_len(n_fam)) {
    g_sub <- igraph::induced_subgraph(g, vids = fam_idx[[i]])

    fam_df <- .fam_df_from_ids(fam_idx[[i]])
    fam_seed <- if (is.null(seed)) NULL else as.integer(seed + i)

    if (nrow(fam_df) <= 5) {
      out <- select_upto5(
        fam_df = fam_df,
        g_sub = g_sub,
        seed = fam_seed,
        return_diagnostics = diagnostics
      )
      selector_used <- "upto5"
    } else {
      out <- select_6plus(
        fam_df = fam_df,
        g_sub = g_sub,
        seed = fam_seed,
        return_diagnostics = diagnostics
      )
      selector_used <- "6plus"
    }

    kept_list[[i]] <- out$kept
    discarded_list[[i]] <- out$discarded

    if (diagnostics) {
      if (!is.null(out$mis_method)) mis_method_vec[i] <- out$mis_method
      uniform_vec[i] <- isTRUE(out$uniform)
      if (!is.null(out$n_max_sets)) n_max_sets_vec[i] <- out$n_max_sets

      node_diag_list[[i]] <- out$node_diagnostics
      fs_selector[i] <- if (is.null(out$selector)) selector_used else out$selector
      if (!is.null(out$structure_name)) fs_structure[i] <- out$structure_name
      fs_n_kept[i] <- length(out$kept)
      fs_n_discarded[i] <- length(out$discarded)
      if (!is.null(out$kept_apriori)) fs_n_kept_apriori[i] <- length(out$kept_apriori)
      if (!is.null(out$discarded_apriori)) fs_n_discarded_apriori[i] <- length(out$discarded_apriori)
      if (!is.null(out$candidate_simple)) fs_n_candidate_simple[i] <- length(out$candidate_simple)
      if (!is.null(out$candidate_complex)) fs_n_candidate_complex[i] <- length(out$candidate_complex)
    }
  }

  kept_all <- unique(unlist(kept_list, use.names = FALSE))
  discarded_all <- unique(unlist(discarded_list, use.names = FALSE))

  node_diagnostics <- NULL
  family_summary <- NULL
  if (diagnostics) {
    # Tables are assembled once here: building them family by family is the main
    # cost of the diagnostics.
    node_diagnostics <- .rkn_rbind_rows(node_diag_list)
    nd_rows <- integer(n_fam)
    for (i in seq_len(n_fam)) {
      if (!is.null(node_diag_list[[i]])) nd_rows[i] <- nrow(node_diag_list[[i]])
    }
    node_diagnostics$family_id <- rep(seq_len(n_fam), nd_rows)
    node_diagnostics$family_size <- rep(unname(fam_sizes), nd_rows)

    family_summary <- .rkn_df(
      family_id = seq_len(n_fam),
      family_size = unname(fam_sizes),
      selector = fs_selector,
      structure_name = fs_structure,
      mis_method = mis_method_vec,
      uniform = uniform_vec,
      n_max_sets = n_max_sets_vec,
      n_kept = fs_n_kept,
      n_discarded = fs_n_discarded,
      n_kept_apriori = fs_n_kept_apriori,
      n_discarded_apriori = fs_n_discarded_apriori,
      n_candidate_simple = fs_n_candidate_simple,
      n_candidate_complex = fs_n_candidate_complex
    )
  }

  .new_rkn_selection(
    kept = kept_all,
    discarded = discarded_all,
    family_sizes = unname(fam_sizes),
    method = method,
    criteria = criteria,
    seed = seed,
    node_diagnostics = node_diagnostics,
    family_summary = family_summary,
    graph = g,
    call = call,
    uniform = uniform
  )
}
