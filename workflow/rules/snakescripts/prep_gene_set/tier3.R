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

df_msmr_strict <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |> 
  filter(probeID %in% gene_universe & phenotype %in% phenotypes) |>
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |> 
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
         .by = "phenotype") |>
  setDT(key = c("biosample", "phenotype", "probeID"))

# coloc
df_coloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.coloc.parquet.gz") %>%
  filter(pheno %in% phenotypes, gene %in% gene_universe) |> 
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric))

# multivariant coloc
df_mvcoloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.mvcoloc.parquet.gz") |> 
  filter(pheno %in% phenotypes, gene %in% gene_universe) |> 
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) |> 
  # take maximum pp h4
  group_by(gene, pheno, biosample) |> 
  slice_max(PP.H4.abf, n = 1, with_ties = FALSE) |> 
  setDT() |> 
  mutate(sig = PP.H4.abf >= 0.8)

# annotate strict mr results
df_msmr_tenk10k[df_msmr_strict,
   `:=`(p_SMR_multi.strict = i.p_SMR_multi,
        lfdr_msmr_pheno.strict = i.lfdr_msmr_pheno),
   on = c("biosample", "phenotype", "probeID")]

# annotate coloc  & mv coloc results
df_msmr_tenk10k[df_coloc, coloc_pph4 := i.PP.H4.abf, on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]
df_msmr_tenk10k[df_mvcoloc, mvcoloc_pph4 := i.PP.H4.abf, on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]

# filter to only significant gene & add tiered evidence
evidence_tier <- list(
  tier1 = expression(lfdr_msmr_pheno < 0.05 & p_HEIDI >= 0.05),
  tier2 = expression(tier1 == TRUE & lfdr_msmr_pheno.strict < 0.05),
  tier3 = expression(tier2 == TRUE & (coloc_pph4 >= 0.8 | mvcoloc_pph4 >= 0.8))
)

# Create columns for each tier
for (tier_name in names(evidence_tier)) {
  df_msmr_tenk10k[, (tier_name) := eval(evidence_tier[[tier_name]])]
}

df_out <- df_msmr_tenk10k[tier3 == TRUE]

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
