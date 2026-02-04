# workflow uso per famiglie >=6

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(here)
})

#rkn_select_6plus_random (rkn_select_6plus_random + apply_rkn_select_6plus_random)
#generate_network_6plus (create_family_graph + make_synthetic_network)
#usf_plots_6plus (plot_full_network + plot_family_diagnostics + plot_full_network_diagnostics + save_families_diagnostics_pdf_grid)

devtools::load_all()
source(here("scripts/generate_network_6plus.R"))
source(here("scripts/usf_plots_6plus.R"))


synth <- make_synthetic_network(n6=25, n7=25, n8=25, n9=25, n10=25, seed=123)
diag  <- apply_rkn_select_random_6plus(synth, seed=999)
plot_full_network(synth, layout = "fr", vertex_size = 3, vertex_label = FALSE, show_legend = F)
plot_full_network_diagnostics(synth, diag, layout="fr", vertex_size=3, vertex_label=FALSE)

plot_family_diagnostics(synth, diag, family_id = 14, layout = "circle", vertex_label = F)

out_pdf <- save_families_diagnostics_pdf_grid(
  synth, diag,
  pdf_path = "rkn_families_diagnostics_6plus.pdf",
  nrow = 8, ncol = 8,
  layout = "fr",
  vertex_size = 15,
  vertex_label = FALSE,
  width = 16, height = 16,
  add_page_title_in_cell = FALSE
)
