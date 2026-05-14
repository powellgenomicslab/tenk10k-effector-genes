
# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(readxl)

OUTPUT <- snakemake@output
df_trait_map <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx") |>
  filter(include)

phenotypes <- df_trait_map$trait_id

# df_msmr_tenk10k <- read_parquet("results/aggregate/tenk10k_phase1.msmr.parquet.gz") %>%
#   filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
#   setDT(key = c("biosample", "phenotype", "probeID"))

df_msmr_tenk10k <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |>
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
  filter(phenotype %in% phenotypes)

# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
gene_universe <- unique(df_msmr_tenk10k$probeID)

writeLines(gene_universe, OUTPUT$gene_universe)
