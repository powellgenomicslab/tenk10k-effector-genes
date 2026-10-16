# Purpose: Count the Kong et al. 2023 DEG genes in the major cell types that carry
# Crohn's disease MR genes, for the "DEG in Matched Major Cell Type" bar of Fig 6a.
# The DEG table itself is not shared, so the bundle carries this count instead,
# together with the set of major cell types it was computed for, which
# scripts/8-crohns/plot_figures/1-gene_annotation_barplot.R checks before using it.

library(tidyverse)

# major cell types as defined in scripts/8-crohns/plot_figures/1-gene_annotation_barplot.R
tenk_sig <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>%
  ungroup() %>%
  mutate(major_cell_type = as.character(major_cell_type) %>% dplyr::replace_when(cell_type == "ILC" ~ "NK", cell_type == "Plasmablast" ~ "Plasma B"))
major_cell_types <- sort(unique(tenk_sig$major_cell_type))

deg <- readRDS("resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS")
deg <- deg %>% filter(major_cell_type %in% major_cell_types)

tibble(
  major_cell_types = paste(major_cell_types, collapse = ";"),
  n_genes = length(unique(deg$Gene))
) %>%
  write_tsv("resources/crohns_case_study/deg/deg_genes_matched_major_ct.tsv")
