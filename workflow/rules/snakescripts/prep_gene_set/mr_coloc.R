# preprocessing script with stringent MSMR results as the baseline results

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(readxl)
library(qvalue)

dir_out <- snakemake@output[[1]]
df_cell_map <- fread("resources/metadata/cell_map.tsv")
df_trait_map_all <- fread("metadata/trait.tsv")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

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
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
         .by = "phenotype") |>
  setDT(key = c("biosample", "phenotype", "probeID"))


# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
df_magma_all <- read_parquet("results/aggregate/tenk10k_phase1.magma.gz.parquet")
# gene_universe <- intersect(df_magma_all$GENE, df_msmr_tenk10k$probeID)
gene_universe <- unique(df_msmr_tenk10k$probeID)

df_msmr_eqtlgen <- read_parquet("results/aggregate/eqtlgen2020.msmr.gz.parquet") %>%
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0,
         probeID %in% gene_universe) %>%
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
         .by = "phenotype") %>%
  mutate(sig = lfdr_msmr_pheno < 0.05 & p_HEIDI >= 0.05) %>%
  setDT()

df_magma <- df_magma_all %>%
  filter(GENE %in% gene_universe, phenotype %in% phenotypes) %>%
  mutate(qval = qvalue(P)$qvalues,
         lfdr = qvalue(P)$lfdr,
         pbh = p.adjust(P, "BH")) %>%
  filter(lfdr < 0.05) %>%
  setDT()

df_msmr_tenk10k <- df_msmr_tenk10k %>%
  # filter(probeID %in% gene_universe) %>%
  # mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
  #        lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
  #        pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
  #        .by = "phenotype") %>%
  # mutate(sig = ifelse(is.na(p_HEIDI), qval_msmr_pheno < 0.05,
  #   qval_msmr_pheno < 0.05 & p_HEIDI >= 0.01
  # )) %>%
  left_join(df_cell_map %>% select(biosample = wg2_scpred_prediction, cell_type, major_cell_type)) %>%
  inner_join(df_trait_map %>%
    select(
      phenotype = trait_id, pheno_label = label,
      pheno_cat = cat_rev, supercategory
    )) %>%
  group_by(phenotype) %>%
  mutate(
    cell_type = factor(cell_type, df_cell_map$cell_type),
    major_cell_type = factor(major_cell_type, unique(df_cell_map$major_cell_type)),
    pheno_cat = factor(pheno_cat, cat_order)
  ) %>%
  setDT(key = c("biosample", "phenotype", "probeID"))

# df_msmr_strict <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |>
#   filter(probeID %in% gene_universe) |>
#   filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
#   mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
#          lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
#          pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
#          .by = "phenotype") |>
#   setDT(key = c("biosample", "phenotype", "probeID"))

# coloc
df_coloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.coloc.parquet.gz") %>%
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) %>%
  mutate(pp_h3_h4 = PP.H3.abf + PP.H4.abf,
         sig = PP.H4.abf >= 0.8) %>%
  filter(pheno %in% phenotypes, gene %in% gene_universe)

# multivariant coloc
df_mvcoloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.mvcoloc.parquet.gz") |>
  filter(pheno %in% phenotypes, gene %in% gene_universe) |>
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) |>
  mutate(pp_h3_h4 = PP.H3.abf + PP.H4.abf) |>
  # take maximum
  group_by(gene, pheno, biosample) |>
  slice_max(PP.H4.abf, n = 1, with_ties = FALSE) |>
  setDT() |>
  mutate(sig = PP.H4.abf >= 0.8)

# annotate gene and phenotype
df_msmr_tenk10k[df_gene_annot, gene_type := i.gene_type, on = c("probeID" = "ensembl_gene_id")]

# annotate strict mr results
# df_msmr_tenk10k[df_msmr_strict,
#    `:=`(p_SMR_multi.strict = i.p_SMR_multi,
#         lfdr_msmr_pheno.strict = i.lfdr_msmr_pheno),
#    on = c("biosample", "phenotype", "probeID")]

# annotate magma results
df_msmr_tenk10k[, magma_gene := FALSE]
df_msmr_tenk10k[df_magma, magma_gene := TRUE, on = c("probeID" = "GENE", "phenotype")]

# annotate eqtlgen mr results
df_msmr_tenk10k[, eqtlgen_mr := FALSE]
df_msmr_tenk10k[df_msmr_eqtlgen, eqtlgen_mr := i.sig, on = c("probeID", "phenotype")]

# annotate coloc  & mv coloc results
df_msmr_tenk10k[df_coloc, coloc_pph4 := i.PP.H4.abf, on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]
df_msmr_tenk10k[df_mvcoloc, mvcoloc_pph4 := i.PP.H4.abf, on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]

evidence_criteria <- list(
  mr = expression(lfdr_msmr_pheno < 0.05 & p_HEIDI >= 0.05),
  coloc = expression(coloc_pph4 >= 0.8 | mvcoloc_pph4 >= 0.8),
  mr_coloc = expression(mr & coloc)
)

# Create columns for each tier
# for (tier_name in names(evidence_tier)) {
#   df_msmr_tenk10k[, (tier_name) := eval(evidence_tier[[tier_name]])]
# }

# Create columns for eviedence criteria and count
for (e in names(evidence_criteria)) {
  df_msmr_tenk10k[, (e) := eval(evidence_criteria[[e]])]
}

df_out <- df_msmr_tenk10k[mr_coloc == TRUE]

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
