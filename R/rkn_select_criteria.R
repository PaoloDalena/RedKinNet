#' Criteria-based selection over an entire network (connected-component-wise)
#'
#' @description
#' Criteria-based counterpart of \code{\link{rkn_select_random}}. The network is
#' partitioned into families, defined as the connected components of \code{g}, and
#' each family is resolved independently: families with at most 5 nodes through the
#' structural rules in \code{\link{rules_db}}, larger families through structural
#' reductions followed by the resolution of the remaining ambiguity. Ambiguous
#' selections are resolved by criterion score rather than uniformly at random; family
#' construction, aggregation of results, console summary, plotting, and return
#' structure are identical to \code{rkn_select_random()}.
#'
#' Covariates referenced by the criteria are read from the vertex attributes of
#' \code{g}.
#'
#' @inheritParams rkn_select_random
#' @param criterion_primary An \code{rkn_criterion} object (see
#'   \code{\link{rkn_criterion_target_mean}} and siblings) used to score candidate
#'   resolutions within each family.
#' @param criterion_secondary Optional \code{rkn_criterion} object used to break ties in
#'   \code{criterion_primary}.
#' @param seed Optional integer base seed. If provided, family-specific seeds are derived
#'   as seed + family_id and are used only to break exact ties in the criterion score(s).
#'
#' @return Same structure as \code{\link{rkn_select_random}}.
#'
#' @seealso \code{\link{rkn_select_random}}, \code{\link{rkn_criterion_target_mean}},
#'   \code{\link{rules_db}}, \code{\link{rkn_plot_pre}}, \code{\link{plot.rkn_selection}}.
#'
#' @export
rkn_select_criteria <- function(g,
                                criterion_primary,
                                criterion_secondary = NULL,
                                seed = NULL,
                                diagnostics = FALSE,
                                return_results_by_family = FALSE,
                                verbose = TRUE,
                                plot = c("no", "pre", "post", "diagnostics"),
                                plot_args = list()) {

  plot <- match.arg(plot)

  if (!igraph::is_igraph(g)) stop("g must be an igraph object.")
  if (igraph::vcount(g) == 0) stop("g must contain at least one vertex.")
  if (is.null(igraph::V(g)$name)) stop("V(g)$name must be defined (vertex identifiers are required).")
  if (!is.list(plot_args)) stop("plot_args must be a named list.")
  if (missing(criterion_primary) || !inherits(criterion_primary, "rkn_criterion")) {
    stop("`criterion_primary` must be an rkn_criterion object (see rkn_criterion_*()).")
  }
  if (!is.null(criterion_secondary) && !inherits(criterion_secondary, "rkn_criterion")) {
    stop("`criterion_secondary` must be an rkn_criterion object or NULL.")
  }

  # Diagnostics plotting needs the diagnostics output; otherwise fall back to the post plot.
  plot_effective <- plot
  if (identical(plot, "diagnostics") && !isTRUE(diagnostics)) {
    warning(
      "plot='diagnostics' requires diagnostics=TRUE. Falling back to plot='post' (kept vs discarded).",
      call. = FALSE
    )
    plot_effective <- "post"
  }

  n_nodes <- igraph::vcount(g)

  # Families are extracted by integer vertex id: extraction by name scans the
  # whole graph for every family and makes the run time quadratic.
  fam_idx <- split(seq_len(n_nodes), igraph::components(g)$membership)
  n_fam <- length(fam_idx)
  fam_sizes <- lengths(fam_idx)

  if (identical(plot_effective, "pre")) {
    do.call(rkn_plot_pre, c(list(g = g), plot_args))
  }

  .pct <- function(x, denom) if (denom == 0) NA_real_ else 100 * x / denom
  .fmt_pct <- function(x, digits = 1) if (is.na(x)) "n/a" else sprintf(paste0("%.", digits, "f%%"), x)
  .fmt_num <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)

  # Extract a per-family data.frame of vertex attributes (robust to missing attributes).
  # Node attributes and the structural measures required by rules_db are computed once
  # on the whole graph: repeating these igraph calls for every family dominates the run
  # time on large networks. Degree and eccentricity (max_path_length) are within-family
  # quantities, since igraph ignores vertices in other components.
  node_cols <- c(list(.name = igraph::V(g)$name), igraph::vertex_attr(g))
  if (is.null(node_cols[["name"]])) node_cols[["name"]] <- node_cols[[".name"]]
  if (is.null(node_cols[["degree"]])) {
    node_cols[["degree"]] <- as.integer(igraph::degree(g, mode = "all"))
  }
  if (is.null(node_cols[["max_path_length"]])) {
    node_cols[["max_path_length"]] <- as.integer(igraph::eccentricity(g, mode = "all"))
  }

  .fam_df_from_ids <- function(ids) {
    cols <- node_cols
    for (j in seq_along(cols)) cols[[j]] <- cols[[j]][ids]
    do.call(.rkn_df, cols)
  }

  if (isTRUE(diagnostics) && isTRUE(return_results_by_family)) {
    results_by_family <- vector("list", n_fam)
  } else {
    results_by_family <- NULL
  }

  kept_list <- vector("list", n_fam)
  discarded_list <- vector("list", n_fam)

  used_largest_ivs_vec <- rep(NA, n_fam)
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
      out <- rkn_select_criteria_upto5(
        fam_df = fam_df,
        criterion_primary = criterion_primary,
        criterion_secondary = criterion_secondary,
        seed = fam_seed,
        return_diagnostics = diagnostics,
        g_sub = g_sub
      )
      selector_used <- "upto5"
    } else {
      out <- rkn_select_criteria_6plus(
        fam_df = fam_df,
        g_sub = g_sub,
        criterion_primary = criterion_primary,
        criterion_secondary = criterion_secondary,
        seed = fam_seed,
        return_diagnostics = diagnostics
      )
      selector_used <- "6plus"
    }

    if (!is.null(results_by_family)) results_by_family[[i]] <- out

    kept_list[[i]] <- out$kept
    discarded_list[[i]] <- out$discarded

    if (diagnostics) {
      used_largest_ivs_vec[i] <- if (!is.null(out$used_largest_ivs)) isTRUE(out$used_largest_ivs) else NA

      node_diag_list[[i]] <- out$node_diagnostics
      fs_selector[i] <- selector_used
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

  if (isTRUE(verbose)) {
    n_upto5 <- sum(fam_sizes <= 5)
    n_6plus <- sum(fam_sizes >= 6)
    nodes_selected <- length(kept_all)

    msg <- c(
      "=== rkn_select_criteria() summary ===",
      sprintf("Total families: %d | size (min/mean/max): %d / %s / %d",
              n_fam, min(fam_sizes), .fmt_num(mean(fam_sizes), 2), max(fam_sizes)),
      sprintf("Families <=5: %d (%s) | Families 6+: %d (%s)",
              n_upto5, .fmt_pct(.pct(n_upto5, n_fam)),
              n_6plus, .fmt_pct(.pct(n_6plus, n_fam))),
      sprintf("Total nodes: %d | Selected nodes: %d (%s)",
              n_nodes, nodes_selected, .fmt_pct(.pct(nodes_selected, n_nodes)))
    )

    if (isTRUE(diagnostics)) {
      n_livs <- sum(used_largest_ivs_vec %in% TRUE, na.rm = TRUE)
      msg <- c(msg, sprintf("Families resolved via largest_ivs: %d (%s)",
                            n_livs, .fmt_pct(.pct(n_livs, n_fam))))
    } else {
      msg <- c(msg, "Families resolved via largest_ivs: n/a (available only with diagnostics=TRUE)")
    }

    msg <- c(msg, "================================")
    cat(paste0(msg, collapse = "\n"), "\n")
  }

  if (identical(plot_effective, "post")) {
    do.call(.rkn_plot_post, c(list(g = g, kept = kept_all), plot_args))
  }

  if (!diagnostics) {
    return(kept_all)
  }

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
    used_largest_ivs = used_largest_ivs_vec,
    n_kept = fs_n_kept,
    n_discarded = fs_n_discarded,
    n_kept_apriori = fs_n_kept_apriori,
    n_discarded_apriori = fs_n_discarded_apriori,
    n_candidate_simple = fs_n_candidate_simple,
    n_candidate_complex = fs_n_candidate_complex
  )

  out_main <- list(
    kept = kept_all,
    discarded = discarded_all,
    node_diagnostics = node_diagnostics,
    family_summary = family_summary
  )

  if (isTRUE(return_results_by_family)) {
    out_main$results_by_family <- results_by_family
  }

  if (identical(plot_effective, "diagnostics")) {
    do.call(.rkn_plot_diagnostics, c(list(g = g, res = out_main), plot_args))
  }

  out_main
}
