# Purpose: Extract the parts of two published supplementary tables that the
# Crohn's case study uses, so that the extracts (not the original files) can be
# shared in the data bundle with a citation.
#
# 1. Kong et al. 2023, Immunity 56:444, doi:10.1016/j.immuni.2023.01.002,
#    single-cell DEG table (sheet 1): immune-cell rows and the discrete-model
#    columns read by scripts/8-crohns/deg/1-prepare_crohns_deg_pre-processed.R,
#    with the Ensembl ID of each gene from the study's feature list
#    (features_original.tsv) joined on. Written as Parquet so values are exact.
# 2. Liu et al. 2023, Nat Genet, doi:10.1038/s41588-023-01384-0,
#    Supplementary Table 8: the distinct gene symbols read by
#    scripts/8-crohns/prepare_data/1-prepare_external_source_annot_boolean_df.R.
#
# Run by rule crohns_third_party_extracts in workflow/rules/misc.smk.

library(tidyverse)
library(arrow)

crohns_dir <- "resources/crohns_case_study"

# Kong et al. 2023 ---------------------------------------------------------------
deg <- readxl::read_xlsx(paste0(crohns_dir, "/deg/Kongetal2023_supplementary/Kongetal2023_Crohns-sc_DEG.xlsx"), sheet = 1)
deg <- deg[, 1:10]
deg <- deg %>%
  filter(grepl(pattern = "Immune", `Cell subset`)) %>%
  select(-`Continuous DE coefficients`, -`Continuous DE coefficients p value`, -`Continuous FDR`)
features <- read_tsv(paste0(crohns_dir, "/deg/features_original.tsv"), col_names = c("probeID", "Gene"))
deg <- deg %>% left_join(features, by = "Gene") %>% select(Gene, probeID, everything())
write_parquet(deg, paste0(crohns_dir, "/deg/kong2023_immune_deg.parquet"))

# Liu et al. 2023 ----------------------------------------------------------------
liu <- readxl::read_xlsx(paste0(crohns_dir, "/crohns_relevant_gene_lists/sources/41588_2023_1384_MOESM4_ESM_Liuetal2023_SupTable.xlsx"), sheet = "ST8", skip = 1) %>%
  select(Gene) %>%
  distinct()
write_tsv(liu, paste0(crohns_dir, "/crohns_relevant_gene_lists/sources/liu2023_st8_genes.tsv"))
