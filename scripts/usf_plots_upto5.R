# ============================================================
# Upto5 isoforms diagnostics PDF grid (style aligned to 6plus)
# - last cell = legend
# - no empty leading cell / no shifting
# - uses rkn_select_random_upto5(return_diagnostics=TRUE)
# ============================================================

library(here)
load("graphs_all.rdata")
devtools::load_all()

suppressPackageStartupMessages({
  library(igraph)
})

# ---- your helper (kept as you wrote) ----
make_fake_df_for_graph <- function(g, comp_id = 1L, iso_id = NA_integer_) {
  deg <- degree(g)
  if (vcount(g) >= 5) {
    dist_mat <- distances(g)
    max_pl <- apply(dist_mat, 1, function(x) max(x[is.finite(x)]))
  } else {
    max_pl <- rep(NA_real_, vcount(g))
  }
  data.frame(
    name = as.character(seq_len(vcount(g))),
    comp_id = comp_id,
    comp_size = vcount(g),
    degree = as.integer(deg),
    iso_id = iso_id,
    max_path_length = as.numeric(max_pl),
    stringsAsFactors = FALSE
  )
}

# ---- color palette: analogous to minimal 6plus palette, but for upto5 classes ----
.rkn_diag_colors_upto5 <- function() {
  c(
    kept_apriori      = "#1B9E77",
    discarded_apriori = "#B3B3B3",
    candidate_simple  = "#FFD92F",
    candidate_complex = "#E41A1C"
  )
}

# ---- compute vertex classes for ONE isoform, from res$node_diagnostics ----
# Priority (strong -> weak): kept_apriori, discarded_apriori, candidate_complex, candidate_simple
# We assign from weaker to stronger so stronger overwrites.
.rkn_vertex_classes_iso_upto5 <- function(g, res) {
  vnames <- V(g)$name
  cls <- setNames(rep("other", length(vnames)), vnames)

  nd <- res$node_diagnostics
  if (is.null(nd) || nrow(nd) == 0) return(cls)

  # expected pre_class values from the updated function
  kept_ap <- intersect(nd$id[nd$pre_class == "kept_apriori"], vnames)
  disc_ap <- intersect(nd$id[nd$pre_class == "discarded_apriori"], vnames)
  csimple <- intersect(nd$id[nd$pre_class == "candidate_simple"], vnames)
  ccomplex <- intersect(nd$id[nd$pre_class == "candidate_complex"], vnames)

  # assign weak -> strong
  if (length(csimple))  cls[csimple]  <- "candidate_simple"
  if (length(ccomplex)) cls[ccomplex] <- "candidate_complex"
  if (length(disc_ap))  cls[disc_ap]  <- "discarded_apriori"
  if (length(kept_ap))  cls[kept_ap]  <- "kept_apriori"

  cls
}

# ---- thumbnail plot (like plot_family_thumbnail_minimal in usf_plots_6plus) ----
plot_isoform_thumbnail_upto5 <- function(g,
                                         res,
                                         layout = c("circle", "fr", "kk"),
                                         vertex_size = 6,
                                         vertex_label = FALSE,
                                         main = NULL) {
  layout <- match.arg(layout)

  cls <- .rkn_vertex_classes_iso_upto5(g, res)
  col_map <- .rkn_diag_colors_upto5()

  vcol <- unname(col_map[cls])
  vcol[is.na(vcol)] <- "#D9D9D9"  # other

  lay <- switch(
    layout,
    circle = layout_in_circle(g),
    fr     = layout_with_fr(g),
    kk     = layout_with_kk(g)
  )

  if (is.null(main)) main <- ""

  plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) V(g)$name else NA,
    vertex.label.cex = 0.5,
    vertex.label.color = "black",
    edge.arrow.size = 0,
    edge.curved = 0.1,
    main = main
  )

  invisible(NULL)
}

# ---- main: save a PDF grid (same logic as save_families_diagnostics_pdf_grid) ----
save_isoforms_diagnostics_pdf_grid_upto5 <- function(graphs_all,
                                                     rules_db,
                                                     pdf_path = "rkn_isoforms_diagnostics_upto5_3to5.pdf",
                                                     iso_names = NULL,
                                                     nrow = 8,
                                                     ncol = 8,
                                                     layout = c("circle", "fr", "kk"),
                                                     vertex_size = 6,
                                                     vertex_label = FALSE,
                                                     width = 16,
                                                     height = 16,
                                                     seed = 123,
                                                     add_page_title_in_cell = TRUE) {
  layout <- match.arg(layout)

  if (is.null(iso_names)) {
    # keep only n3..n5 isoforms if names follow "n3_isoX" pattern
    nm_all <- names(graphs_all)
    iso_names <- nm_all[grepl("^n[3-5]_iso", nm_all)]
  }

  N <- length(iso_names)
  if (N == 0) stop("No isoforms found. Check graphs_all names or pass iso_names explicitly.")

  total_panels <- nrow * ncol
  plots_per_page <- total_panels - 1  # last = legend
  n_pages <- ceiling(N / plots_per_page)

  col_map <- .rkn_diag_colors_upto5()

  pdf(pdf_path, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  idx_global <- 0L

  for (p in seq_len(n_pages)) {
    start_idx <- (p - 1) * plots_per_page + 1
    end_idx <- min(p * plots_per_page, N)
    these <- iso_names[start_idx:end_idx]

    par(
      mfrow = c(nrow, ncol),
      mar = c(0.2, 0.2, 0.8, 0.2),
      oma = c(0.3, 0.3, 0.3, 0.3)
    )

    # Fill exactly plots_per_page cells with isoforms (or blank), then legend cell.
    for (i in seq_len(plots_per_page)) {
      if (i <= length(these)) {
        nm <- these[i]
        g <- graphs_all[[nm]]

        # iso_id from name like "n4_iso26"
        iso_id <- suppressWarnings(as.integer(sub(".*iso", "", nm)))
        if (is.na(iso_id)) iso_id <- NA_integer_

        fam_df <- make_fake_df_for_graph(g, comp_id = 1L, iso_id = iso_id)

        res <- rkn_select_random_upto5(
          fam_df = fam_df,
          rules_db = rules_db,
          seed = seed,
          return_diagnostics = TRUE,
          g_sub = g
        )

        star <- if (isTRUE(res$used_largest_ivs)) " *" else ""
        if (isTRUE(add_page_title_in_cell) && i == 1) {
          main <- paste0("page ", p, "/", n_pages, " - ", nm, " -> ", res$structure_name, star)
        } else {
          main <- paste0(nm, " -> ", res$structure_name, star)
        }

        plot_isoform_thumbnail_upto5(
          g = g,
          res = res,
          layout = layout,
          vertex_size = vertex_size,
          vertex_label = vertex_label,
          main = main
        )
      } else {
        plot.new()
      }
    }

    # Last cell = legend (same idea as 6plus)
    plot.new()
    legend(
      "center",
      legend = names(col_map),
      col = unname(col_map),
      pch = 19,
      pt.cex = 1.2,
      cex = 1.5,
      bty = "n",
      title = "Diagnostic class:"
    )
  }

  message("Saved PDF: ", pdf_path)
  invisible(pdf_path)
}

# -------------------
# Example usage:
# -------------------
save_isoforms_diagnostics_pdf_grid_upto5(
  graphs_all = graphs_all,
  rules_db = rules_db,
  pdf_path = "rkn_upto5_isoforms_3to5.pdf",
  nrow = 6, ncol = 5,
  layout = "fr",
  vertex_size = 15,
  seed = 123
)
