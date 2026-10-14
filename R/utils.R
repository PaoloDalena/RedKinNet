#' Internal helper: lightweight data frame constructor
#'
#' Builds a data frame from equal-length atomic vectors without the argument
#' processing of \code{data.frame()}, whose overhead dominates when thousands of
#' small tables are created (one or more per family). Character columns stay
#' character and column names are dropped from the values.
#'
#' @keywords internal
#' @noRd
.rkn_df <- function(...) {
  x <- list(...)
  for (j in seq_along(x)) names(x[[j]]) <- NULL
  n <- if (length(x) == 0) 0L else length(x[[1]])
  structure(x, class = "data.frame", row.names = .set_row_names(n))
}

#' Internal helper: row-bind a list of data frames with identical columns
#'
#' Concatenates column by column, matching columns by name on the first data frame;
#' NULL entries are ignored.
#'
#' @keywords internal
#' @noRd
.rkn_rbind_rows <- function(dfs) {
  keep <- logical(length(dfs))
  for (k in seq_along(dfs)) keep[k] <- !is.null(dfs[[k]])
  dfs <- dfs[keep]
  if (length(dfs) == 0) return(NULL)
  cols <- names(dfs[[1]])
  out <- vector("list", length(cols))
  names(out) <- cols
  for (j in seq_along(cols)) {
    pieces <- vector("list", length(dfs))
    for (k in seq_along(dfs)) pieces[[k]] <- .subset2(dfs[[k]], cols[j])
    out[[j]] <- unlist(pieces, use.names = FALSE)
  }
  do.call(.rkn_df, out)
}
