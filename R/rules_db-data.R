#' Rule database for connected components up to 5 nodes
#'
#' A curated list of structural rules used by \code{\link{rkn_select_random}} and
#' \code{\link{rkn_select_criteria}} on families with at most 5 nodes, to perform
#' deterministic (a priori) keep/discard decisions and to define conflict groups for
#' the resolution of the remaining ambiguity. Every connected graph with at most 5
#' nodes matches exactly one rule.
#'
#' Each element of the list is a rule with fields such as:
#' \itemize{
#'   \item \code{name}
#'   \item \code{size_allowed}
#'   \item \code{match_fn}
#'   \item \code{kept_rules}
#'   \item \code{discarded_rules}
#'   \item \code{conflict_rules}
#'   \item optional \code{use_largest_ivs} (logical): if \code{TRUE}, the remaining
#'     ambiguity is resolved by choosing among all maximum independent sets of the
#'     eligible nodes
#'   \item optional \code{kept_fn}, \code{discarded_fn}
#' }
#'
#' @format A list of rule objects (lists).
#'
#' @examples
#' length(rules_db)
#' rules_db[[2]]$name
"rules_db"
