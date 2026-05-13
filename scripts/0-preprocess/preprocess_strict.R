# preprocessing script with stringent MSMR results as the baseline results

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(readxl)
library(qvalue)

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

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v3.parquet.gz") 

# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
df_magma_all <- read_parquet("results/aggregate/tenk10k_phase1.magma.gz.parquet")

gene_universe <- unique(df_msmr_tenk10k$probeID)

df_msmr_eqtlgen <- read_parquet("results/aggregate/eqtlgen2020.msmr.parquet.gz") %>%
  # impute p_SMR_multi if missing but p_SMR is available
  mutate(p_SMR_multi = ifelse(is.na(p_SMR_multi) & !is.na(p_SMR), p_SMR, p_SMR_multi)) %>%
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0,
         probeID %in% gene_universe) %>% 
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
         .by = "phenotype") %>% 
  mutate(sig = lfdr_msmr_pheno < 0.05) %>% 
  setDT()

df_magma <- df_magma_all %>%
  filter(GENE %in% gene_universe, phenotype %in% phenotypes) %>%
  mutate(qval = qvalue(P)$qvalues,
         lfdr = qvalue(P)$lfdr,
         pbh = p.adjust(P, "BH")) %>%
  filter(lfdr < 0.05) %>%
  setDT()

# annotate magma results
df_msmr_tenk10k[, magma_gene := FALSE]
df_msmr_tenk10k[df_magma, magma_gene := TRUE, on = c("probeID" = "GENE", "phenotype")]

# annotate eqtlgen mr results
df_msmr_tenk10k[, eqtlgen_mr := FALSE]
df_msmr_tenk10k[df_msmr_eqtlgen, eqtlgen_mr := i.sig, on = c("probeID", "phenotype")]

