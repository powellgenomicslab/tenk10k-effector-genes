source("scripts/util/write_table.R")
library(tidyverse)

# Read spearman correlation results (written by scripts/2-trait/fig_2_heatmap_spearman_calculation.R)
# Columns: p1, p2, spearman_corr, n_intersecting_genes
df_spearman <- read_tsv("results/rrho/spearman_corr_all_trait_combos_strictmr.tsv")

# column labels come from the spearman_corr sheet of metadata/table_column_names.xlsx:
#     ~name, ~label,
#     "p1", "Phenotype 1",
#     "p2", "Phenotype 2",
#     "spearman_corr", "Spearman Correlation",
#     "n_intersecting_genes", "Number of Intersecting Genes"

write_table(df_spearman, source_tbl_name = "spearman_corr", tbl_order = 4)
