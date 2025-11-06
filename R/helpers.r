# Helper: apply .node_filter safely
apply_rule_filter <- function(fam_df, rule_exprs) {
  if (length(rule_exprs) == 0) return(character(0))
  out <- character(0)
  for (expr_txt in rule_exprs) {
    mask <- .node_filter(fam_df, expr_txt)
    # se per qualsiasi ragione arriva uno scalare, replicalo
    if (length(mask) == 1L) {
      mask <- rep(isTRUE(mask), nrow(fam_df))
    } else {
      mask <- as.logical(mask)
      mask[is.na(mask)] <- FALSE
    }
    idx <- which(mask)
    out <- c(out, fam_df$.name[idx])
  }
  unique(out)
}

# Helper: find which rule matches a given component
find_matching_rule <- function(fam_df, rules_db) {
  for (rule in rules_db) {
    if (rule$size_allowed == nrow(fam_df) && isTRUE(rule$match_fn(fam_df))) {
      return(rule)
    }
  }
  return(NULL)
}

.node_filter <- function(fam_df, expr_txt) {
  # gestisci i casi banali in modo vettoriale
  if (identical(expr_txt, "TRUE"))  return(rep(TRUE,  nrow(fam_df)))
  if (identical(expr_txt, "FALSE")) return(rep(FALSE, nrow(fam_df)))
  # usa un env "lista" per vedere le colonne come vettori
  eval(parse(text = expr_txt), envir = as.list(fam_df))
}
