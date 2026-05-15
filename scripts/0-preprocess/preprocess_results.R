# preprocess main results

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(arrow)
  library(fs)
  library(readxl)
  library(qvalue)
})

major_cell_type_order <- c(
  "CD4 T", "CD8 T", "Unconventional T", "NK", "Plasma B", "B", "Monocyte", "Dendritic", "HSPC"
)

df_cell_map <- fread("resources/metadata/cell_map_revised.tsv") |>
  select(-major_cell_type) |>
  rename(major_cell_type = revision_major_cell_type) |>
  mutate(major_cell_type = factor(major_cell_type, levels = major_cell_type_order))

df_trait_map_all <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

setDT(df_trait_map)
phenotypes <- df_trait_map$trait_id

cat_order <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx",
  sheet = "trait_category_order"
) %>%
  pull(cat_order)

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v4.parquet.gz")

df_magma_all <- read_parquet("results/aggregate/tenk10k_phase1.magma.gz.parquet")

gene_universe <- unique(df_msmr_tenk10k$probeID)

df_msmr_eqtlgen <- read_parquet("results/aggregate/eqtlgen2020.msmr.parquet.gz") %>%
  mutate(p_SMR_multi = ifelse(is.na(p_SMR_multi) & !is.na(p_SMR), p_SMR, p_SMR_multi)) %>%
  filter(
    !is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0,
    probeID %in% gene_universe
  ) %>%
  mutate(
    qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
    lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
    pbh_msmr_pheno  = p.adjust(p_SMR_multi, "BH"),
    .by = "phenotype"
  ) %>%
  mutate(sig = lfdr_msmr_pheno < 0.05) %>%
  setDT()

df_magma <- df_magma_all %>%
  filter(GENE %in% gene_universe, phenotype %in% phenotypes) %>%
  mutate(
    qval = qvalue(P)$qvalues,
    lfdr = qvalue(P)$lfdr,
    pbh  = p.adjust(P, "BH")
  ) %>%
  filter(lfdr < 0.05) %>%
  setDT()

# annotate MAGMA and eQTLGen overlap onto TenK10K results
df_msmr_tenk10k[, magma_gene := FALSE]
df_msmr_tenk10k[df_magma, magma_gene := TRUE, on = c("probeID" = "GENE", "phenotype")]

df_msmr_tenk10k[, eqtlgen_mr := FALSE]
df_msmr_tenk10k[df_msmr_eqtlgen, eqtlgen_mr := i.sig, on = c("probeID", "phenotype")]

# coloc
df_coloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.coloc.parquet.gz") %>%
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) %>%
  mutate(
    pp_h3_h4 = PP.H3.abf + PP.H4.abf,
    sig = PP.H4.abf >= 0.8
  ) %>%
  filter(pheno %in% phenotypes, gene %in% gene_universe)

# multivariant coloc
df_mvcoloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.mvcoloc.parquet.gz") |>
  filter(pheno %in% phenotypes, gene %in% gene_universe) |>
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) |>
  mutate(pp_h3_h4 = PP.H3.abf + PP.H4.abf) |>
  group_by(gene, pheno, biosample) |>
  slice_max(PP.H4.abf, n = 1, with_ties = FALSE) |>
  setDT() |>
  mutate(sig = PP.H4.abf >= 0.8)

# significant MR results
df_msmr <- df_msmr_tenk10k[mr == TRUE]
pheno_order <- df_msmr[, .N, by = pheno_label][order(-N), pheno_label]
df_msmr[, pheno_label := factor(pheno_label, pheno_order)]
