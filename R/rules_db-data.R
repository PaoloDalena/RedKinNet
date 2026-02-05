#' Rule database for connected components up to 5 nodes
#'
#' A curated list of structural rules used by \code{\link{rkn_select_random_upto5}}
#' to perform deterministic (a priori) keep/discard decisions and to define
#' conflict groups for randomized resolution.
#'
#' Each element of the list is a rule with fields such as:
#' \itemize{
#'   \item \code{name}
#'   \item \code{size_allowed}
#'   \item \code{match_fn}
#'   \item \code{kept_rules}
#'   \item \code{discarded_rules}
#'   \item \code{conflict_rules}
#'   \item optional \code{alpha}, \code{kept_fn}, \code{discarded_fn}
#' }
#'
#' @format A list of rule objects (lists).
#'
#' @examples
#' # Inspect available structure names ruled
#' sapply(rules_db, "[[", 1)
#'
"rules_db"
