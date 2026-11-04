#' Internal helper: criteria-based selection over all maximum sets of one family
#'
#' Chooses the best maximum set of unrelated individuals of the family by the criteria,
#' with \code{.rkn_criteria_maximum_set()} applied to the whole family instead of the
#' individuals left undecided by the structural rules, which fix some individuals
#' regardless of the criteria. With diagnostics, the exact inclusion probabilities over
#' the maximum sets searched classify each individual as in every maximum set
#' (\code{"kept_apriori"}), in no maximum set (\code{"discarded_apriori"}) or decided by
#' the criteria among alternative sets (\code{"candidate_complex"}). A family too
#' intricate for the exact engine, even after the dominance reduction, is selected with
#' the structural rules instead, as in the default selection, and reported with
#' \code{optimal = FALSE}.
#'
#' @return The per-family result read by \code{.rkn_select_network()}, with
#'   \code{optimal} and, with diagnostics, \code{n_max_sets} (\code{NA} after the
#'   dominance reduction, which leaves out some maximum sets).
#'
#' @keywords internal
#' @noRd
.rkn_select_optimal_family <- function(fam_df, g_sub, criterion_primary, criterion_secondary, seed,
                                       return_diagnostics, max_width = 20) {
  ids <- igraph::V(g_sub)$name
  picked <- .rkn_criteria_maximum_set(g_sub, fam_df, character(0), criterion_primary, criterion_secondary,
                                      seed = seed, max_width = max_width, error_if_wide = FALSE)

  # Families of up to 5 members always fit within the width limit.
  if (is.null(picked)) {
    out <- rkn_select_criteria_6plus(fam_df = fam_df, g_sub = g_sub, criterion_primary = criterion_primary,
                                     criterion_secondary = criterion_secondary, seed = seed,
                                     return_diagnostics = return_diagnostics)
    out$optimal <- FALSE
    return(out)
  }

  kept <- picked$kept
  discarded <- setdiff(ids, kept)
  if (!return_diagnostics) return(list(kept = kept, discarded = discarded, optimal = TRUE))

  searched <- .rkn_search_inclusion(g_sub, picked$searched)
  p <- searched$inclusion
  every <- ids[p == 1]
  none <- ids[p == 0]
  chosen <- ids[p > 0 & p < 1]
  reason <- stats::setNames(rep("optimal_chosen", length(ids)), ids)
  reason[every] <- "optimal_every"
  reason[none] <- "optimal_none"
  node_diag <- .rkn_build_node_diag_6plus(
    all_names = ids, kept_apriori = every, discarded_apriori = none,
    candidate_simple = character(0), candidate_complex = chosen,
    kept_final = kept, discarded_final = discarded, pre_reason_map = reason
  )

  list(
    kept = kept, discarded = discarded, optimal = TRUE, selector = "optimal",
    structure_name = NA_character_, mis_method = picked$method,
    kept_apriori = every, discarded_apriori = none,
    candidate_simple = character(0), candidate_complex = chosen,
    n_max_sets = if (identical(picked$method, "exact_reduced")) NA_real_ else searched$n_max_sets,
    node_diagnostics = node_diag
  )
}
