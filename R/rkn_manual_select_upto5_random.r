rkn_manual_select_upto5_random <- function(fam_df,
                                           rules_db,
                                           seed = NULL,
                                           return_diagnostics = FALSE,
                                           g_sub = NULL) {
  # 0) allinea name / .name
  if (is.null(fam_df$.name) && !is.null(fam_df$name)) fam_df$.name <- fam_df$name
  if (is.null(fam_df$name)  && !is.null(fam_df$.name)) fam_df$name  <- fam_df$.name

  # 1) trova la regola
  rule <- find_matching_rule(fam_df, rules_db)

  # helper: costruisci diagnostica nodo-per-nodo
  build_node_diag <- function(ids,
                              kept_apriori, discarded_apriori,
                              candidate_simple, candidate_complex,
                              kept_final, discarded_final,
                              pre_reason_map = NULL) {
    pre_class <- rep(NA_character_, length(ids))
    names(pre_class) <- ids

    pre_class[ids %in% kept_apriori]      <- "kept_apriori"
    pre_class[ids %in% discarded_apriori] <- "discarded_apriori"
    pre_class[ids %in% candidate_simple]  <- "candidate_simple"
    pre_class[ids %in% candidate_complex] <- "candidate_complex"

    # fallback: tutto ciò che resta (non dovrebbe succedere) -> candidate_complex
    pre_class[is.na(pre_class)] <- "candidate_complex"

    post_class <- rep(NA_character_, length(ids))
    names(post_class) <- ids
    # assegna per nome (robusto anche se kept_final/discarded_final sono id in character)
    post_class[intersect(ids, kept_final)]      <- "kept"
    post_class[intersect(ids, discarded_final)] <- "discarded"
    post_class[is.na(post_class)] <- "discarded"

    pre_reason <- rep(NA_character_, length(ids))
    names(pre_reason) <- ids
    if (!is.null(pre_reason_map) && length(pre_reason_map) > 0) {
      pre_reason[names(pre_reason_map)] <- unname(pre_reason_map)
    }

    data.frame(
      id = ids,
      pre_class = unname(pre_class[ids]),
      post_class = unname(post_class[ids]),
      is_selected = unname(post_class[ids] == "kept"),
      pre_reason = unname(pre_reason[ids]),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(rule)) {
    kept_final <- character(0)
    discarded_final <- fam_df$.name

    if (!return_diagnostics) {
      return(list(kept = kept_final, discarded = discarded_final))
    }

    node_diagnostics <- data.frame(
      id = fam_df$.name,
      pre_class = rep("discarded_apriori", nrow(fam_df)),
      post_class = rep("discarded", nrow(fam_df)),
      is_selected = rep(FALSE, nrow(fam_df)),
      pre_reason = rep("no_matching_rule", nrow(fam_df)),
      stringsAsFactors = FALSE
    )

    family_diagnostics <- data.frame(
      n_nodes = nrow(fam_df),
      used_largest_ivs = FALSE,
      structure_name = "unknown",
      stringsAsFactors = FALSE
    )

    return(list(
      structure_name = "unknown",
      kept_apriori = character(0),
      discarded_apriori = character(0),
      candidate_simple = character(0),
      candidate_complex = character(0),
      kept_afterrandomization = character(0),
      discarded_afterrandomization = character(0),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = FALSE,
      node_diagnostics = node_diagnostics,
      family_diagnostics = family_diagnostics
    ))
  }

  # 2) deterministici tabellari
  kept_apriori      <- apply_rule_filter(fam_df, unlist(rule$kept_rules))
  discarded_apriori <- apply_rule_filter(fam_df, unlist(rule$discarded_rules))

  # 2b) deterministici grafo-dipendenti (kept_fn / discarded_fn)
  if (!is.null(g_sub)) {
    if (!is.null(rule$kept_fn) && is.function(rule$kept_fn)) {
      kept_graph <- rule$kept_fn(fam_df, g_sub)
      kept_apriori <- unique(c(kept_apriori, kept_graph))
    }
    if (!is.null(rule$discarded_fn) && is.function(rule$discarded_fn)) {
      disc_graph <- rule$discarded_fn(fam_df, g_sub)
      discarded_apriori <- unique(c(discarded_apriori, disc_graph))
    }
  }

  # 3) gruppi di conflitto (escludendo SEMPRE gli a-priori)
  conflict_groups <- list()
  candidates_all  <- character(0)
  if (length(rule$conflict_rules) > 0) {
    for (grp in rule$conflict_rules) {
      cand_nodes <- apply_rule_filter(fam_df, grp)
      cand_nodes <- setdiff(cand_nodes, union(kept_apriori, discarded_apriori))
      cand_nodes <- unique(cand_nodes)
      if (length(cand_nodes) > 0) conflict_groups[[length(conflict_groups)+1]] <- cand_nodes
    }
    if (length(conflict_groups) > 0) {
      candidates_all <- unique(unlist(conflict_groups))
    }
  }

  # 4) selezione: largest_ivs se alpha>1 e ho g_sub, altrimenti random classico
  used_largest_ivs <- FALSE
  kept_afterrandomization <- character(0)
  discarded_afterrandomization <- character(0)

  if (!is.null(rule$alpha) && isTRUE(rule$alpha > 1) && !is.null(g_sub)) {
    # eleggibili = tutti i nodi tranne discarded a priori e (opzionale) vicini dei kept a priori
    eligible <- setdiff(fam_df$.name, discarded_apriori)

    if (length(kept_apriori) > 0) {
      nb_idx  <- igraph::neighborhood(g_sub, order = 1, nodes = kept_apriori)
      nb_names <- unique(setdiff(unlist(lapply(nb_idx, function(vs) igraph::V(g_sub)$name[vs])), kept_apriori))
      eligible <- setdiff(eligible, union(kept_apriori, nb_names))
    }

    g_ind <- igraph::induced_subgraph(g_sub, vids = eligible)
    sets  <- igraph::largest_ivs(g_ind)
    set_names <- lapply(sets, function(s) igraph::V(g_ind)$name[s])
    if (!is.null(seed)) set.seed(seed)
    chosen_names <- if (length(set_names) > 0) set_names[[sample(seq_along(set_names), 1)]] else character(0)

    kept_afterrandomization <- unique(chosen_names)
    # se non ho definito candidates con i gruppi, usa eligible per le diagnostiche
    if (length(candidates_all) == 0) candidates_all <- eligible
    discarded_afterrandomization <- setdiff(candidates_all, kept_afterrandomization)
    used_largest_ivs <- TRUE

  } else if (length(conflict_groups) > 0) {
    if (!is.null(seed)) set.seed(seed)
    if (length(conflict_groups) == 1) {
      grp <- conflict_groups[[1]]
      kept_afterrandomization <- if (length(grp) == 1) grp else sample(grp, size = 1)
    } else {
      U <- unique(unlist(conflict_groups))
      k <- length(conflict_groups)
      if (length(U) >= k) {
        kept_afterrandomization <- sample(U, size = k, replace = FALSE)
      } else {
        chosen <- character(0)
        for (grp in conflict_groups) {
          avail <- setdiff(grp, chosen)
          if (length(avail) > 0) chosen <- c(chosen, sample(avail, 1))
        }
        kept_afterrandomization <- unique(chosen)
      }
    }
    discarded_afterrandomization <- setdiff(candidates_all, kept_afterrandomization)
  }

  kept_afterrandomization <- unique(kept_afterrandomization)

  # 5) finali
  kept_final <- unique(c(kept_apriori, kept_afterrandomization))
  discarded_final <- unique(c(discarded_apriori, discarded_afterrandomization))
  discarded_final <- setdiff(discarded_final, kept_final)  # kept vince

  # 6) complessità (candidate_complex se più gruppi OPPURE alpha>1)
  is_complex <- (length(conflict_groups) > 1) ||
    (!is.null(rule$alpha) && isTRUE(rule$alpha > 1))

  candidate_complex <- if (is_complex) unique(candidates_all) else character(0)
  candidate_simple  <- if (!is_complex) unique(candidates_all) else character(0)

  if (!return_diagnostics) {
    return(list(kept = kept_final, discarded = discarded_final))
  } else {
    # pre_reason minimale (opzionale ma utile)
    pre_reason_map <- character(0)
    if (used_largest_ivs && length(candidate_complex) > 0) {
      pre_reason_map[candidate_complex] <- "largest_ivs_fallback"
    } else if (length(conflict_groups) > 1 && length(candidate_complex) > 0) {
      pre_reason_map[candidate_complex] <- "multi_conflict_groups"
    } else if (length(conflict_groups) == 1 && length(candidate_simple) > 0) {
      pre_reason_map[candidate_simple] <- "single_conflict_group"
    }

    node_diagnostics <- build_node_diag(
      ids = fam_df$.name,
      kept_apriori = unique(kept_apriori),
      discarded_apriori = unique(discarded_apriori),
      candidate_simple = unique(candidate_simple),
      candidate_complex = unique(candidate_complex),
      kept_final = kept_final,
      discarded_final = discarded_final,
      pre_reason_map = pre_reason_map
    )

    family_diagnostics <- data.frame(
      n_nodes = nrow(fam_df),
      used_largest_ivs = used_largest_ivs,
      structure_name = rule$name,
      stringsAsFactors = FALSE
    )

    return(list(
      structure_name = rule$name,
      kept_apriori = unique(kept_apriori),
      discarded_apriori = unique(discarded_apriori),
      candidate_simple = unique(candidate_simple),
      candidate_complex = unique(candidate_complex),
      kept_afterrandomization = unique(kept_afterrandomization),
      discarded_afterrandomization = unique(discarded_afterrandomization),
      kept = kept_final,
      discarded = discarded_final,
      used_largest_ivs = used_largest_ivs,
      node_diagnostics = node_diagnostics,
      family_diagnostics = family_diagnostics
    ))
  }
}
