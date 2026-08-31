# preprocessing script with stringent MSMR results as the baseline results

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
# library(qvalue)
library(readxl)

dir_out <- snakemake@output[[1]]
df_cell_map <- fread("resources/metadata/cell_map.tsv")
df_trait_map_all <- fread("metadata/trait.tsv")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

# STRING ID mapping (ensembl_gene_id → stringId)
df_string_map <- fread(snakemake@input$string_id_map)
cat("STRING ID map:", nrow(df_string_map), "entries loaded\n")

setDT(df_trait_map)
phenotypes <- df_trait_map$trait_id

cat_order <- fread("metadata/trait_category.tsv") %>%
  pull(cat_order)

# df_msmr_tenk10k <- read_parquet("results/aggregate/tenk10k_phase1.msmr.parquet.gz") %>%
#   filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
#   setDT(key = c("biosample", "phenotype", "probeID"))

df_msmr_tenk10k <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |>
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
  filter(phenotype %in% phenotypes) |>
  # mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
  #        lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
  #        pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
  #        .by = "phenotype") |>
  left_join(df_cell_map %>% select(biosample = wg2_scpred_prediction, cell_type, major_cell_type)) %>%
  inner_join(df_trait_map %>%
    select(
      phenotype = trait_id, pheno_label = label,
      pheno_cat = cat_rev, supercategory
    )) %>%
  mutate(
    cell_type = factor(cell_type, df_cell_map$cell_type),
    major_cell_type = factor(major_cell_type, unique(df_cell_map$major_cell_type)),
    pheno_cat = factor(pheno_cat, cat_order)
  ) %>%
  setDT(key = c("biosample", "phenotype", "probeID"))

# annotate gene and phenotype
df_msmr_tenk10k[df_gene_annot, gene_type := i.gene_type, on = c("probeID" = "ensembl_gene_id")]

df_msmr_tenk10k[df_string_map, stringId := i.stringId, on = c("probeID" = "ensembl_gene_id")]

for (b in c(unique(df_msmr_tenk10k$biosample), "overall")) {
  if (b == "overall") {
    df <- df_msmr_tenk10k
  } else {
    df <- df_msmr_tenk10k[biosample == b]
  }
  if (nrow(df) == 0) {
    next
  }
  for (p in unique(df$phenotype)) {
    # Signed statistic: sign(b_SMR) * -log10(p_SMR_multi)
    # Positive = upregulated & significant, Negative = downregulated & significant
    # This enables STRING's both-direction enrichment detection
    df_out <- df[phenotype == p & !is.na(stringId), list(
      stringId,
      score = sign(b_SMR) * (-log10(p_SMR_multi))
    )][order(-abs(score))]
    if (nrow(df_out) == 0) {
      next
    }
    outfile <- file.path(dir_out, b, paste0(p, ".txt"))
    dir_create(dirname(outfile))
    fwrite(df_out, outfile, sep = "\t", col.names = FALSE)
  }
}
