#' Internal helper: uniform selection within one family
#'
#' Draws one maximum set of unrelated individuals of the family uniformly at random
#' among all of them, with the exact engine applied to the whole family (no structural
#' rule and no dominance reduction, which would make some maximum sets unreachable).
#' With diagnostics, the exact inclusion probability of each individual classifies it as
#' in every maximum set (\code{"kept_apriori"}), in no maximum set
#' (\code{"discarded_apriori"}) or chosen at random among maximum sets
#' (\code{"candidate_complex"}). A family whose elimination width exceeds
#' \code{max_width} is selected by the structural rules instead, as in the default
#' selection, and reported with \code{uniform = FALSE}.
#'
#' @return The per-family result read by \code{.rkn_select_network()}, with
#'   \code{uniform} and, with diagnostics, \code{n_max_sets} and the column
#'   \code{inclusion_probability}.
#'
#' @keywords internal
#' @noRd
.rkn_select_uniform_family <- function(fam_df, g_sub, seed, return_diagnostics, max_width = 20) {
  ids <- igraph::V(g_sub)$name
  sol <- .rkn_mis_exact(.rkn_adjacency(g_sub), max_width)

  if (!sol$feasible) {
    out <- if (length(ids) <= 5) {
      rkn_select_random_upto5(fam_df = fam_df, seed = seed, return_diagnostics = return_diagnostics, g_sub = g_sub)
    } else {
      rkn_select_random_6plus(fam_df = fam_df, g_sub = g_sub, seed = seed, return_diagnostics = return_diagnostics)
    }
    out$uniform <- FALSE
    if (return_diagnostics) out$node_diagnostics$inclusion_probability <- NA_real_
    return(out)
  }

  if (!is.null(seed)) {
    withr::local_preserve_seed()
    set.seed(seed)
  }
  x <- .rkn_mis_sample(sol)
  kept <- ids[x]
  discarded <- ids[!x]
  if (!return_diagnostics) return(list(kept = kept, discarded = discarded, uniform = TRUE))

  p <- .rkn_mis_inclusion(sol)
  every <- ids[p == 1]
  none <- ids[p == 0]
  chosen <- ids[p > 0 & p < 1]
  reason <- stats::setNames(rep("uniform_chosen", length(ids)), ids)
  reason[every] <- "uniform_every"
  reason[none] <- "uniform_none"
  node_diag <- .rkn_build_node_diag_6plus(
    all_names = ids, kept_apriori = every, discarded_apriori = none,
    candidate_simple = character(0), candidate_complex = chosen,
    kept_final = kept, discarded_final = discarded, pre_reason_map = reason
  )
  node_diag$inclusion_probability <- p

  list(
    kept = kept, discarded = discarded, uniform = TRUE, selector = "uniform",
    structure_name = NA_character_, mis_method = "exact",
    kept_apriori = every, discarded_apriori = none,
    candidate_simple = character(0), candidate_complex = chosen,
    n_max_sets = round(exp(sol$log_count)), node_diagnostics = node_diag
  )
}
