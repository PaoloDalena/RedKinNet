library(here)
source(here("R", "helpers.r"))
source(here("R", "rkn_manual_select_upto5_random.r"))
source(here("R", "rules_db.r"))
load(here("graphs_all.rdata"))

make_fake_df_for_graph <- function(g, comp_id = 1L, iso_id = NA_integer_) {
  deg <- degree(g)
  if (vcount(g) >= 5) {
    dist_mat <- distances(g)
    max_pl <- apply(dist_mat, 1, function(x) max(x[is.finite(x)]))
  } else {
    max_pl <- rep(NA_real_, vcount(g))
  }
  data.frame(name = as.character(seq_len(vcount(g))), comp_id = comp_id, comp_size = vcount(g), degree = as.integer(deg), iso_id = iso_id, max_path_length = as.numeric(max_pl), stringsAsFactors = FALSE)
}

plot_graph_basic <- function(g, main = NULL) {
  V(g)$label <- V(g)$name
  V(g)$size <- 26
  V(g)$frame.color <- "black"
  V(g)$label.color <- "black"
  E(g)$color <- "grey40"
  E(g)$width <- 1.2
  lay <- if (vcount(g) > 2 && ecount(g) == vcount(g) - 1 && max(degree(g)) <= 2) {
    layout_as_tree(g, circular = FALSE)
  } else if (vcount(g) == 1) {
    matrix(c(0, 0), ncol = 2)
  } else {
    layout_in_circle(g)
  }
  plot(g, layout = lay, main = main)
}

plot_graph_colored <- function(g, colors, main = NULL) {
  V(g)$color <- colors[V(g)$name]
  plot_graph_basic(g, main = main)
}

# --- loop sui grafi (aggiornato) ---
sel <- names(graphs_all)

old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

pdf("diagnostica_redkinnet_upto5.pdf", width = 12, height = 6)

for (nm in sel) {
  g <- graphs_all[[nm]]

  # estraggo iso_id dal nome tipo "n4_iso26"
  iso_id <- suppressWarnings(as.integer(sub(".*iso", "", nm)))
  if (is.na(iso_id)) iso_id <- NA_integer_

  df_fake <- make_fake_df_for_graph(g, comp_id = 1L, iso_id = iso_id)

  cat("\n===========================\n")
  cat("Grafico:", nm, "\n")
  print(df_fake)

  par(mfrow = c(1, 2), mar = c(2,2,3,1), oma = c(0,0,2,0))

  plot_graph_basic(g, main = paste0(nm, " (raw)"))

  res <- rkn_manual_select_upto5_random(
    fam_df = df_fake,
    rules_db = rules_db,
    seed = 123,
    return_diagnostics = TRUE,
    g_sub = g         # <— NOVITÀ: serve per largest_ivs()
  )

  # colori
  node_names <- df_fake$name
  cols <- setNames(rep("white", length(node_names)), node_names)

  if (length(res$candidates) > 0) {
    if (length(res$candidates_complex) > 0) {
      cols[ res$candidates ] <- "orangered"
    } else {
      cols[ res$candidates ] <- "yellow"
    }
  }
  if (length(res$discarded_apriori) > 0) {
    cols[ res$discarded_apriori ] <- "grey80"
  }
  if (length(res$kept_apriori) > 0) {
    cols[ res$kept_apriori ] <- "forestgreen"
  }

  # titolo con asterisco se scelto via largest_ivs
  star <- if (isTRUE(res$used_largest_ivs)) " *" else ""
  plot_graph_colored(
    g,
    colors = cols,
    main = paste0(nm, " \u2192 ", res$structure_name, star)
  )
}

dev.off()
