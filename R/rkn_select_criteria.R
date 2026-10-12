#' Criteria-based selection over an entire network (connected-component-wise)
#'
#' @description
#' Criteria-based counterpart of \code{\link{rkn_select_random}}. The network is
#' partitioned into families, defined as the connected components of \code{g}, and
#' each family is resolved independently by \code{\link{rkn_select_criteria_upto5}} (families with at
#' most 5 nodes) or \code{\link{rkn_select_criteria_6plus}} (families with 6 or more
#' nodes). Ambiguous selections are resolved by criterion score rather than uniformly at
#' random; family construction, aggregation of results, console summary, plotting, and
#' return structure are identical to \code{rkn_select_random()}.
#'
#' Covariates referenced by the criteria are read from the vertex attributes of
#' \code{g}, which are propagated to the per-family data frames. The packaged
#' \code{rules_db} is always used for families with at most 5 nodes.
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
#' @seealso \code{\link{rkn_select_random}}, \code{\link{rkn_select_criteria_upto5}},
#'   \code{\link{rkn_select_criteria_6plus}}, \code{rkn_plot_pre()}, \code{rkn_plot_post()},
#'   \code{rkn_plot_diagnostics()}.
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

  all_vnames <- igraph::V(g)$name
  n_nodes <- igraph::vcount(g)

  fam_nodes <- split(all_vnames, igraph::components(g)$membership)
  n_fam <- length(fam_nodes)
  fam_sizes <- lengths(fam_nodes)

  if (identical(plot_effective, "pre")) {
    do.call(rkn_plot_pre, c(list(g = g), plot_args))
  }

  .pct <- function(x, denom) if (denom == 0) NA_real_ else 100 * x / denom
  .fmt_pct <- function(x, digits = 1) if (is.na(x)) "n/a" else sprintf(paste0("%.", digits, "f%%"), x)
  .fmt_num <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)

  # Extract a per-family data.frame of vertex attributes (robust to missing attributes).
  .fam_df_from_graph <- function(g_sub) {
    vn <- igraph::V(g_sub)$name
    attrs <- igraph::vertex_attr(g_sub)
    df <- data.frame(.name = vn, stringsAsFactors = FALSE)
    if (length(attrs) > 0) {
      for (nm in names(attrs)) df[[nm]] <- attrs[[nm]]
    }
    if (is.null(df$name)) df$name <- df$.name

    # Ensure structural measures required by rules_db are present
    if (is.null(df$degree)) {
      df$degree <- as.integer(igraph::degree(g_sub, mode = "all"))
    }
    if (is.null(df$max_path_length)) {
      # Node eccentricity = maximum shortest-path length from each node within the component
      df$max_path_length <- as.integer(igraph::eccentricity(g_sub, mode = "all"))
    }
    df
  }

  if (isTRUE(diagnostics) && isTRUE(return_results_by_family)) {
    results_by_family <- vector("list", n_fam)
  } else {
    results_by_family <- NULL
  }

  kept_all <- character(0)
  discarded_all <- character(0)

  used_largest_ivs_vec <- rep(NA, n_fam)
  node_diag_list <- if (diagnostics) vector("list", n_fam) else NULL
  family_summary_list <- if (diagnostics) vector("list", n_fam) else NULL

  for (i in seq_len(n_fam)) {
    nodes_i <- fam_nodes[[i]]

    g_sub <- igraph::induced_subgraph(g, vids = nodes_i)
    if (is.null(igraph::V(g_sub)$name)) igraph::V(g_sub)$name <- nodes_i

    fam_df <- .fam_df_from_graph(g_sub)
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

    kept_all <- unique(c(kept_all, out$kept))
    discarded_all <- unique(c(discarded_all, out$discarded))

    if (diagnostics) {
      used_largest_ivs_vec[i] <- if (!is.null(out$used_largest_ivs)) isTRUE(out$used_largest_ivs) else NA

      nd <- out$node_diagnostics
      if (!is.null(nd) && nrow(nd) > 0) {
        nd$family_id <- i
        nd$family_size <- nrow(fam_df)
      }
      node_diag_list[[i]] <- nd

      family_summary_list[[i]] <- data.frame(
        family_id = i,
        family_size = nrow(fam_df),
        selector = selector_used,
        structure_name = if (!is.null(out$structure_name)) out$structure_name else NA_character_,
        used_largest_ivs = used_largest_ivs_vec[i],
        n_kept = length(out$kept),
        n_discarded = length(out$discarded),
        n_kept_apriori = if (!is.null(out$kept_apriori)) length(out$kept_apriori) else NA_integer_,
        n_discarded_apriori = if (!is.null(out$discarded_apriori)) length(out$discarded_apriori) else NA_integer_,
        n_candidate_simple = if (!is.null(out$candidate_simple)) length(out$candidate_simple) else NA_integer_,
        n_candidate_complex = if (!is.null(out$candidate_complex)) length(out$candidate_complex) else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

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
    do.call(rkn_plot_post, c(list(g = g, kept = kept_all), plot_args))
  }

  if (!diagnostics) {
    return(kept_all)
  }

  out_main <- list(
    kept = kept_all,
    discarded = discarded_all,
    node_diagnostics = .rkn_rbind_rows(node_diag_list),
    family_summary = .rkn_rbind_rows(family_summary_list)
  )

  if (isTRUE(return_results_by_family)) {
    out_main$results_by_family <- results_by_family
  }

  if (identical(plot_effective, "diagnostics")) {
    do.call(rkn_plot_diagnostics, c(list(g = g, res = out_main), plot_args))
  }

  out_main
}
