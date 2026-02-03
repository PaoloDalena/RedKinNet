#' Internal helper: apply a set of rule expressions to a component data frame
#'
#' @description
#' Evaluates one or more logical expressions (as text) in the context of
#' \code{fam_df} and returns the corresponding node names (from \code{fam_df$.name})
#' where the expression evaluates to \code{TRUE}.
#'
#' Expressions are evaluated via \code{.node_filter()}, which treats the columns
#' of \code{fam_df} as vectors in the evaluation environment.
#'
#' @param fam_df A data frame describing a connected component. Must contain a
#'   character column \code{.name}.
#' @param rule_exprs Character vector of expressions. Each element must evaluate
#'   to a logical vector of length \code{nrow(fam_df)}, or to a scalar logical
#'   (which will be recycled).
#'
#' @return A character vector of unique node names matching at least one rule.
#'
#' @details
#' The function enforces a fully vectorized output mask:
#' \itemize{
#'   \item scalar logical results are recycled to \code{nrow(fam_df)};
#'   \item \code{NA} values are coerced to \code{FALSE}.
#' }
#'
#' @keywords internal
#' @noRd
apply_rule_filter <- function(fam_df, rule_exprs) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (is.null(fam_df$.name)) stop("`fam_df` must contain a `.name` column.")
  if (length(rule_exprs) == 0) return(character(0))

  out <- character(0)

  for (expr_txt in rule_exprs) {
    mask <- .node_filter(fam_df, expr_txt)

    if (length(mask) == 1L) {
      mask <- rep(isTRUE(mask), nrow(fam_df))
    } else {
      mask <- as.logical(mask)
      mask[is.na(mask)] <- FALSE
    }

    out <- c(out, fam_df$.name[which(mask)])
  }

  unique(out)
}

#' Internal helper: find the first rule matching a component
#'
#' @description
#' Scans \code{rules_db} and returns the first rule whose size constraint matches
#' \code{nrow(fam_df)} and whose \code{match_fn(fam_df)} returns \code{TRUE}.
#'
#' @param fam_df A data frame describing a connected component.
#' @param rules_db A list of rule objects. Each rule must provide:
#'   \itemize{
#'     \item \code{size_allowed}: integer.
#'     \item \code{match_fn}: function returning a scalar logical.
#'   }
#'
#' @return The first matching rule (a list), or \code{NULL} if no rule matches.
#'
#' @keywords internal
#' @noRd
find_matching_rule <- function(fam_df, rules_db) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (!is.list(rules_db)) stop("`rules_db` must be a list.")

  for (rule in rules_db) {
    if (isTRUE(rule$size_allowed == nrow(fam_df)) && isTRUE(rule$match_fn(fam_df))) {
      return(rule)
    }
  }

  NULL
}

#' Internal helper: evaluate a rule expression in the context of a component
#'
#' @description
#' Evaluates \code{expr_txt} (a string) as an R expression in an environment where
#' the columns of \code{fam_df} are available as vectors. Two trivial cases are
#' treated explicitly to guarantee vectorized outputs.
#'
#' @param fam_df A data frame describing a connected component.
#' @param expr_txt A single character string containing an R expression.
#'
#' @return The raw evaluation result (typically a logical vector). The caller is
#'   responsible for coercion and NA handling.
#'
#' @details
#' This function uses \code{eval(parse())}. It is intended for controlled inputs
#' coming from an internal rule database, not for untrusted user-provided strings.
#'
#' @keywords internal
#' @noRd
.node_filter <- function(fam_df, expr_txt) {
  if (!is.data.frame(fam_df)) stop("`fam_df` must be a data.frame.")
  if (!is.character(expr_txt) || length(expr_txt) != 1L) {
    stop("`expr_txt` must be a single character string.")
  }

  if (identical(expr_txt, "TRUE"))  return(rep(TRUE,  nrow(fam_df)))
  if (identical(expr_txt, "FALSE")) return(rep(FALSE, nrow(fam_df)))

  eval(parse(text = expr_txt), envir = as.list(fam_df))
}
