# prepare gene set for enrichment analysis

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(readxl)
library(qvalue)

dir_out <- snakemake@output[[1]]
setDTthreads(Sys.getenv("NCPUS"))
df_trait_map <- fread("metadata/trait.tsv") |> 
  filter(include)
phenotypes <- unique(df_trait_map$trait_id)

df_msmr_tenk10k <- read_parquet("results/aggregate/tenk10k_phase1.msmr.parquet.gz") %>%
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
  setDT(key = c("biosample", "phenotype", "probeID"))

# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
df_magma_all <- read_parquet("results/aggregate/tenk10k_phase1.magma.gz.parquet")
gene_universe <- intersect(df_magma_all$GENE, df_msmr_tenk10k$probeID)

df_msmr_tenk10k <- df_msmr_tenk10k %>%
  filter(probeID %in% gene_universe & phenotype %in% phenotypes) %>%
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
         .by = "phenotype")


# filter to only significant gene & add tiered evidence
evidence_tier <- list(
  tier1 = expression(lfdr_msmr_pheno < 0.05 & p_HEIDI >= 0.05)
)

# Create columns for each tier
for (tier_name in names(evidence_tier)) {
  df_msmr_tenk10k[, (tier_name) := eval(evidence_tier[[tier_name]])]
}

df_out <- df_msmr_tenk10k[tier1 == TRUE]

for (b in c(unique(df_out$biosample), "overall")) {
  if (b == "overall") {
    df <- df_out
  } else {
    df <- df_out[biosample == b]
  }
  if (nrow(df) == 0) {
    next
  }
  for (p in unique(df$phenotype)) {
    gene_set <- df[phenotype == p, unique(probeID)]
    if (length(gene_set) == 0) {
      next
    }
    outfile <- file.path(dir_out, b, paste0(p, ".txt"))
    dir_create(dirname(outfile))
    writeLines(gene_set, outfile)
  }
}
