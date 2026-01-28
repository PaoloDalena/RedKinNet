# ============================================================
# PLOT UTILITIES (works with your existing: synth + diag objects)
#   - plot_full_network(): plots whole synthetic network
#   - plot_family_diagnostics(): plots 1 family with diagnostic colors
# ============================================================

# ---------- helper: stable color palette for family_ids (full network)
.rkn_palette <- function(n) {
  grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[1:n],
    c = 70,
    l = 55
  )
}

# ---------- 1) Plot the whole synthetic network (all families)
plot_full_network <- function(synth,
                              layout = c("fr", "kk"),
                              vertex_size = 4,
                              vertex_label = FALSE,
                              edge_arrow_size = 0,
                              show_legend = TRUE,
                              main = "Synthetic network (all families)") {
  layout <- match.arg(layout)

  g <- synth$g_all
  if (is.null(g)) stop("synth$g_all not found.")

  fam_ids <- igraph::V(g)$family_id
  if (is.null(fam_ids)) stop("igraph::V(g)$family_id missing. Use make_synthetic_network().")

  fam_levels <- sort(unique(fam_ids))
  pal <- .rkn_palette(length(fam_levels))
  fam_col_map <- stats::setNames(pal, fam_levels)

  vcol <- unname(fam_col_map[as.character(fam_ids)])

  lay <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g)
  )

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  graphics::par(mar = c(0.5, 0.5, 2, 0.5))

  igraph::plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) igraph::V(g)$name else NA,
    vertex.label.cex = 0.6,
    edge.arrow.size = edge_arrow_size,
    edge.curved = 0.1,
    main = main
  )

  if (isTRUE(show_legend)) {
    # keep legend compact: show at most 25 items
    if (length(fam_levels) <= 25) {
      graphics::legend(
        "topleft",
        legend = paste0("F", fam_levels),
        col = pal,
        pch = 19,
        cex = 0.7,
        bty = "n",
        title = "family_id"
      )
    } else {
      graphics::mtext(
        side = 3,
        line = 0.2,
        adj = 1,
        cex = 0.8,
        paste0("Legend hidden (", length(fam_levels), " families). Set show_legend=FALSE or reduce families.")
      )
    }
  }

  invisible(NULL)
}

# ---------- helper: build a vertex class vector from diag for 1 family
.rkn_vertex_classes_for_family <- function(diag, family_id) {
  # diag$results_by_family is a tibble with columns family_id and out (list)
  if (is.null(diag$results_by_family)) stop("diag$results_by_family not found.")
  row_idx <- which(diag$results_by_family$family_id == family_id)
  if (length(row_idx) != 1) stop("family_id not found (or duplicated) in diag$results_by_family.")

  out <- diag$results_by_family$out[[row_idx]]

  # priority matters for color assignment (top wins):
  # kept_apriori > kept > discarded_apriori > discarded > unknown > complex > simple > undecided
  list(
    structure_name = out$structure_name,
    kept_apriori = out$kept_apriori %||% character(0),
    kept = out$kept %||% character(0),
    discarded_apriori = out$discarded_apriori %||% character(0),
    discarded = out$discarded %||% character(0),
    unknown = out$unknown %||% character(0),
    complex_candidates = out$complex_candidates %||% character(0),
    simple_candidates = out$simple_candidates %||% character(0),
    used_largest_ivs = out$used_largest_ivs %||% FALSE
  )
}

# small infix helper (avoid purrr dependency)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------- 2) Plot a single family with diagnostic colors
plot_family_diagnostics <- function(synth,
                                    diag,
                                    family_id,
                                    layout = c("fr", "kk", "circle"),
                                    vertex_size = 18,
                                    vertex_label = TRUE,
                                    show_legend = TRUE,
                                    main = NULL) {
  layout <- match.arg(layout)

  if (is.null(synth$family_objects)) stop("synth$family_objects not found.")
  idx <- which(synth$family_objects$family_id == family_id)
  if (length(idx) != 1) stop("family_id not found (or duplicated) in synth$family_objects.")

  g <- synth$family_objects$g_sub[[idx]]
  if (is.null(igraph::V(g)$name)) stop("Family g_sub has no vertex names.")

  info <- .rkn_vertex_classes_for_family(diag, family_id)
  structure_name <- info$structure_name
  used_mis <- isTRUE(info$used_largest_ivs)

  vnames <- igraph::V(g)$name
  cls <- rep("undecided", length(vnames))
  names(cls) <- vnames

  # apply in increasing priority, then overwrite with higher priority later
  if (length(info$simple_candidates))  cls[info$simple_candidates] <- "simple_candidates"
  if (length(info$complex_candidates)) cls[info$complex_candidates] <- "complex_candidates"
  if (length(info$unknown))            cls[info$unknown] <- "unknown"
  if (length(info$discarded))          cls[info$discarded] <- "discarded"
  if (length(info$discarded_apriori))  cls[info$discarded_apriori] <- "discarded_apriori"
  if (length(info$kept))               cls[info$kept] <- "kept"
  if (length(info$kept_apriori))       cls[info$kept_apriori] <- "kept_apriori"

  # color map (feel free to change)
  col_map <- c(
    kept_apriori = "#1B9E77",
    kept = "#66C2A5",
    discarded_apriori = "#D95F02",
    discarded = "#FDAE61",
    simple_candidates = "#7570B3",
    complex_candidates = "#E7298A",
    unknown = "#666666",
    undecided = "#BDBDBD"
  )

  vcol <- unname(col_map[cls])

  lay <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g),
    circle = igraph::layout_in_circle(g)
  )

  if (is.null(main)) {
    main <- paste0(
      "Family ", family_id,
      " | ", structure_name,
      if (used_mis) " | largest_ivs=TRUE" else " | largest_ivs=FALSE"
    )
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, 2.2, 0.5))

  igraph::plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) vnames else NA,
    vertex.label.cex = 0.75,
    vertex.label.color = "black",
    edge.arrow.size = 0,
    edge.curved = 0.1,
    main = main
  )

  if (isTRUE(show_legend)) {
    present <- intersect(names(col_map), unique(cls))
    graphics::legend(
      "topleft",
      legend = present,
      col = unname(col_map[present]),
      pch = 19,
      cex = 0.8,
      bty = "n",
      title = "diagnostic class"
    )
  }

  invisible(list(class = cls, colors = vcol))
}


.rkn_diag_colors <- function() {
  c(
    kept_apriori = "#1B9E77",
    discarded_apriori = "#B3B3B3",
    simple_candidates = "#FFD92F",
    complex_candidates = "#E41A1C",
    unknown = "#4D4D4D",
    other = "#D9D9D9"
  )
}

# Costruisce una mappa name -> class usando diag$results_by_family$out
# Classi richieste: kept_apriori, discarded_apriori, simple_candidates, complex_candidates, unknown
# (Per robustezza: se un nodo non cade in nessuna classe, va in "other")
rkn_build_node_class_map <- function(synth, diag) {
  if (is.null(synth$g_all)) stop("synth$g_all non trovato.")
  if (is.null(diag$results_by_family)) stop("diag$results_by_family non trovato.")

  all_names <- igraph::V(synth$g_all)$name
  cls_map <- stats::setNames(rep("other", length(all_names)), all_names)

  # helper: assegna classe a un vettore di nodi
  assign_class <- function(nodes, class_name) {
    nodes <- intersect(nodes, all_names)
    if (length(nodes) > 0) cls_map[nodes] <<- class_name
  }

  # Qui definiamo la PRIORITÀ di assegnazione (l’ultima vince).
  # Vogliamo che kept/discarded a priori dominino sugli altri.
  for (i in seq_len(nrow(diag$results_by_family))) {
    out <- diag$results_by_family$out[[i]]

    kept_ap <- out$kept_apriori %||% character(0)
    disc_ap <- out$discarded_apriori %||% character(0)
    simple  <- out$simple_candidates %||% character(0)
    complex <- out$complex_candidates %||% character(0)
    unknown <- out$unknown %||% character(0)

    # ordine: prima metto "più debole", poi sovrascrivo con "più forte"
    assign_class(simple,  "simple_candidates")
    assign_class(complex, "complex_candidates")
    assign_class(unknown, "unknown")
    assign_class(disc_ap, "discarded_apriori")
    assign_class(kept_ap, "kept_apriori")
  }

  cls_map
}

# ------------------------------------------------------------
# 1) Plot full network colorato per diagnostica
# ------------------------------------------------------------
plot_full_network_diagnostics <- function(synth,
                                          diag,
                                          layout = c("fr", "kk"),
                                          vertex_size = 4,
                                          vertex_label = FALSE,
                                          show_legend = TRUE,
                                          main = "Full network (colored by diagnostics)") {
  layout <- match.arg(layout)

  g <- synth$g_all
  if (is.null(g)) stop("synth$g_all non trovato.")
  if (is.null(igraph::V(g)$name)) stop("igraph::V(g)$name mancante in g_all.")

  cls_map <- rkn_build_node_class_map(synth, diag)
  cls <- unname(cls_map[igraph::V(g)$name])

  col_map <- .rkn_diag_colors()
  vcol <- unname(col_map[cls])
  vcol[is.na(vcol)] <- col_map["other"]

  lay <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g)
  )

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, 2, 0.5))

  igraph::plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) igraph::V(g)$name else NA,
    vertex.label.cex = 0.6,
    vertex.label.color = "black",
    edge.arrow.size = 0,
    edge.curved = 0.1,
    main = main
  )

  if (isTRUE(show_legend)) {
    present <- sort(unique(cls))
    present <- intersect(present, names(col_map))
    graphics::legend(
      "topleft",
      legend = present,
      col = unname(col_map[present]),
      pch = 19,
      cex = 0.85,
      bty = "n",
      title = "Diagnostic class:"
    )
  }

  invisible(NULL)
}

# ============================================================
# PDF GRID (FIXED): no cella vuota, legenda sempre nell'ultima cella
# Classi SOLO come in plot_full_network_diagnostics:
#   kept_apriori, discarded_apriori, simple_candidates, complex_candidates, unknown (+ other)
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.rkn_diag_colors_minimal <- function() {
  c(
    kept_apriori = "#1B9E77",
    discarded_apriori = "#B3B3B3",
    simple_candidates = "#FFD92F",
    complex_candidates = "#E41A1C",
    unknown = "#4D4D4D"
  )
}

.rkn_get_out_for_family <- function(diag, family_id) {
  if (is.null(diag$results_by_family)) stop("diag$results_by_family non trovato.")
  idx <- which(diag$results_by_family$family_id == family_id)
  if (length(idx) != 1) stop("family_id non trovato (o duplicato) in diag$results_by_family.")
  diag$results_by_family$out[[idx]]
}

# Classificazione MINIMA (5 classi) coerente con plot_full_network_diagnostics
# Priorità (forte -> debole): kept_apriori, discarded_apriori, unknown, complex, simple
.rkn_vertex_classes_family_minimal <- function(g_family, out) {
  vnames <- igraph::V(g_family)$name
  cls <- stats::setNames(rep("other", length(vnames)), vnames)

  kept_ap <- intersect(out$kept_apriori %||% character(0), vnames)
  disc_ap <- intersect(out$discarded_apriori %||% character(0), vnames)
  simple  <- intersect(out$simple_candidates %||% character(0), vnames)
  complex <- intersect(out$complex_candidates %||% character(0), vnames)
  unknown <- intersect(out$unknown %||% character(0), vnames)

  # assegna dal più debole al più forte
  if (length(simple))  cls[simple]  <- "simple_candidates"
  if (length(complex)) cls[complex] <- "complex_candidates"
  if (length(unknown)) cls[unknown] <- "unknown"
  if (length(disc_ap)) cls[disc_ap] <- "discarded_apriori"
  if (length(kept_ap)) cls[kept_ap] <- "kept_apriori"

  cls
}

# Thumbnail plot: NON modifica par(), NON mette legenda, classi minime
plot_family_thumbnail_minimal <- function(synth,
                                          diag,
                                          family_id,
                                          layout = c("circle", "fr", "kk"),
                                          vertex_size = 6,
                                          vertex_label = FALSE,
                                          main = NULL) {
  layout <- match.arg(layout)

  idx <- which(synth$family_objects$family_id == family_id)
  if (length(idx) != 1) stop("family_id non trovato (o duplicato) in synth$family_objects.")
  g <- synth$family_objects$g_sub[[idx]]
  if (is.null(igraph::V(g)$name)) stop("g_sub senza V$name.")

  out <- .rkn_get_out_for_family(diag, family_id)
  cls <- .rkn_vertex_classes_family_minimal(g, out)

  col_map <- .rkn_diag_colors_minimal()
  vcol <- unname(col_map[cls])
  vcol[is.na(vcol)] <- col_map["other"]

  lay <- switch(
    layout,
    circle = igraph::layout_in_circle(g),
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g)
  )

  if (is.null(main)) main <- paste0("F", family_id)

  igraph::plot.igraph(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size = vertex_size,
    vertex.label = if (isTRUE(vertex_label)) igraph::V(g)$name else NA,
    vertex.label.cex = 0.5,
    vertex.label.color = "black",
    edge.arrow.size = 0,
    edge.curved = 0.1,
    main = main
  )

  invisible(NULL)
}

# PDF grid: ultimo riquadro = legenda, nessun buco.
save_families_diagnostics_pdf_grid <- function(synth,
                                               diag,
                                               pdf_path = "rkn_families_diagnostics_8x8.pdf",
                                               family_ids = NULL,
                                               nrow = 8,
                                               ncol = 8,
                                               layout = c("circle", "fr", "kk"),
                                               vertex_size = 6,
                                               vertex_label = FALSE,
                                               width = 16,
                                               height = 16,
                                               add_page_title_in_cell = TRUE) {
  layout <- match.arg(layout)

  if (is.null(family_ids)) family_ids <- synth$family_objects$family_id

  total_panels <- nrow * ncol
  plots_per_page <- total_panels - 1
  N <- length(family_ids)
  n_pages <- ceiling(N / plots_per_page)

  col_map <- .rkn_diag_colors_minimal()

  grDevices::pdf(file = pdf_path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (p in seq_len(n_pages)) {
    start_idx <- (p - 1) * plots_per_page + 1
    end_idx <- min(p * plots_per_page, N)
    these_ids <- family_ids[start_idx:end_idx]

    graphics::par(
      mfrow = c(nrow, ncol),
      mar = c(0.2, 0.2, 0.8, 0.2),
      oma = c(0.3, 0.3, 0.3, 0.3)
    )

    # Riempio esattamente (plots_per_page) celle con famiglie (o vuote se mancano),
    # poi l'ultima cella è legenda.
    for (i in seq_len(plots_per_page)) {
      if (i <= length(these_ids)) {
        fid <- these_ids[i]
        # Se vuoi il titolo pagina senza rompere mfrow, scrivilo SOLO nella prima cella
        if (isTRUE(add_page_title_in_cell) && i == 1) {
          # titolo piccolo dentro la cella (in alto), senza unicode “—”
          plot_family_thumbnail_minimal(
            synth = synth, diag = diag, family_id = fid,
            layout = layout, vertex_size = vertex_size, vertex_label = vertex_label,
            main = paste0("page ", p, "/", n_pages, " - F", fid)
          )
        } else {
          plot_family_thumbnail_minimal(
            synth = synth, diag = diag, family_id = fid,
            layout = layout, vertex_size = vertex_size, vertex_label = vertex_label,
            main = paste0("F", fid)
          )
        }
      } else {
        graphics::plot.new()
      }
    }

    # Ultima cella: legenda
    graphics::plot.new()
    graphics::legend(
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

  invisible(pdf_path)
}
