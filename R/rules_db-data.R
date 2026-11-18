#' Structural rules for families of up to 5 members
#'
#' A curated list of structural rules used by \code{\link{rkn_select_random}} and
#' \code{\link{rkn_select_criteria}} on families of up to 5 members: each rule
#' recognises a family structure, keeps or discards some members outright, and groups
#' the others into the alternatives left to choose from. Every possible structure of a
#' family of up to 5 members matches exactly one rule. The rules are not used in uniform
#' selections (\code{rkn_select_random(uniform = TRUE)}) nor in selections over all
#' maximum sets (\code{rkn_select_criteria(optimal = TRUE)}).
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
#'     ambiguity is resolved by choosing among all the maximum sets of the eligible
#'     members
#'   \item optional \code{kept_fn}, \code{discarded_fn}
#' }
#'
#' @format A list of rule objects (lists).
#'
#' @examples
#' length(rules_db)
#' rules_db[[2]]$name
"rules_db"
