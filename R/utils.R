#' Internal helper: row-bind a list of data frames, ignoring NULL entries
#'
#' @keywords internal
#' @noRd
.rkn_rbind_rows <- function(dfs) {
  out <- do.call(rbind, dfs)
  if (!is.null(out)) rownames(out) <- NULL
  out
}
