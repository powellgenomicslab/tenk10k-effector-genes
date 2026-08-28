# preprocessing script with stringent MSMR results as the baseline results
# version 5: 2026-08-10
# Update with N instruments used in MR
# Manually replace ASD results
# Update IVW and MRlink2 results with v2 

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(readxl)
library(qvalue)

setDTthreads(Sys.getenv("NCPUS"))
df_cell_map <- fread("resources/metadata/cell_map_revised.tsv") |> 
  select(-major_cell_type) |> 
  rename(major_cell_type = revision_major_cell_type)
df_trait_map_all <- fread("metadata/trait.tsv")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

setDT(df_trait_map)
phenotypes <- df_trait_map$trait_id

cat_order <- fread("metadata/trait_category.tsv") %>%
  pull(cat_order)

df_msmr_tenk10k <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |> 
  filter(phenotype %in% phenotypes) |>
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |> 
  filter(phenotype %in% phenotypes) |> 
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         .by = "phenotype") |> 
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

# annotate gene and phenotype
df_msmr_tenk10k[df_gene_annot, gene_type := i.gene_type, on = c("probeID" = "ensembl_gene_id")]

gene_universe <- unique(df_msmr_tenk10k$probeID)

# mrlink2 - explicit source path (results/aggregate/ holds only a symlink to this)
df_mrlink2 <- read_parquet("results/sensitivity/mrlink2/tenk10k_phase1_sensitivity.mrlink2.v2.parquet.gz")

# ivw-ld - the `robust` variant; results/aggregate/ holds only a symlink to this
df_ivw <- read_parquet("results/sensitivity/ivw-ld/robust/tenk10k_phase1_sensitivity.ivw-ld.v2.parquet.gz")

# coloc - v3: rebuilt 2026-08-13 after fixing the aggregation reader and
# re-running coloc for the 34 traits whose GWAS carried zero-SE or duplicate-ID
# SNPs. See CHANGELOG.v5.md.
df_coloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.coloc.v3.parquet.gz") %>%
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) %>%
  mutate(pp_h3_h4 = PP.H3.abf + PP.H4.abf,
         sig = PP.H4.abf >= 0.8) %>%
  filter(pheno %in% phenotypes, gene %in% gene_universe)

# multivariant coloc
df_mvcoloc <- read_parquet("results/aggregate/coloc/tenk10k_phase1.mvcoloc.parquet.gz") |> 
  filter(pheno %in% phenotypes, gene %in% gene_universe) |> 
  mutate(across(c(PP.H0.abf:PP.H4.abf), as.numeric)) |> 
  mutate(pp_h3_h4 = PP.H3.abf + PP.H4.abf) |> 
  # take maximum PP.H4.abf per gene-phenotype-biosample combination
  group_by(gene, pheno, biosample) |>
  arrange(desc(PP.H4.abf), .by_group = TRUE) |> 
  slice_head(n = 1) |> 
  setDT()

# annotate ivw-ld results
df_msmr_tenk10k[df_ivw, `:=`(
  b_ivw = i.estimate,
  se_ivw = i.se,
  p_ivw = i.pval,
  phet_ivw = i.het_pval,
  fstat_ivw = i.fstat,
  isq_ivw = i.i_sq
  ),
  on = c("biosample", "phenotype", "probeID")]

# annotate mrlink2 results
df_msmr_tenk10k[df_mrlink2, `:=`(
    b_mrlink2 = i.alpha,
    se_mrlink2 = `i.se(alpha)`,
    p_mrlink2 = `i.p(alpha)`,
    psigmay_mrlink2 = `i.p(sigma_y)`,
    sesigmay_mrlink2 = `i.se(sigma_y)`,
    sigmay_mrlink2 = `i.sigma_y` 
  ), on = c("biosample", "phenotype", "probeID")]

# annotate coloc  & mv coloc results
df_msmr_tenk10k[df_coloc, `:=`(coloc_pph3_h4 = i.pp_h3_h4, coloc_pph4 = i.PP.H4.abf), on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]
df_msmr_tenk10k[df_mvcoloc, `:=`(mvcoloc_pph4 = i.PP.H4.abf), on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]

# arrange and rename columns

cols <- c("biosample", "cell_type", "major_cell_type",
          "phenotype", "pheno_label", "pheno_cat", "supercategory",
          "probeID", "Gene", "gene_type", "ProbeChr", "Probe_bp",
          "topSNP", "topSNP_chr", "topSNP_bp",
          "A1", "A2", "Freq",
          "b_GWAS", "se_GWAS", "p_GWAS",
          "b_eQTL", "se_eQTL", "p_eQTL",
          "b_SMR", "se_SMR", "p_SMR", "p_SMR_multi", "p_HEIDI", "nsnp_HEIDI", "n_snps_smr",
          "qval_msmr_pheno", "lfdr_msmr_pheno",
          "b_ivw", "se_ivw", "p_ivw", "phet_ivw", 
          "fstat_ivw", "isq_ivw",
          "b_mrlink2", "se_mrlink2", "p_mrlink2", "sigmay_mrlink2", "sesigmay_mrlink2", "psigmay_mrlink2",
          "coloc_pph3_h4", "coloc_pph4", "mvcoloc_pph4")

df_msmr_tenk10k <- df_msmr_tenk10k[, ..cols]

# NA handling in the OR'd criteria.
#
# R's `|` already resolves TRUE | NA to TRUE, so a single positive source wins
# over missing ones - that part was always correct. The problem is the other
# direction: `FALSE | NA` is NA, so a known-negative source combined with a
# missing one produced NA rather than FALSE. Because mvcoloc_pph4 is NA for
# 93.6% of rows (mvcoloc legitimately finds nothing for most gene x trait
# pairs), 93% of `coloc` values were NA - of which 7,005,035 had a known
# coloc_pph4 < 0.8 and should have been FALSE. Downstream that meant
# `filter(coloc == FALSE)` returned 420,612 rows instead of 7,425,647, since
# NA == FALSE is NA.
#
# Each comparison is therefore wrapped in `%in% TRUE`, which maps NA to FALSE,
# i.e. a missing source is treated as non-supporting. NA is preserved only where
# *no* source was assessed at all, so "not supported" stays distinguishable from
# "not tested" (11,487 rows for coloc).
evidence_criteria <- list(
  mr = expression(lfdr_msmr_pheno < 0.05),
  sensitivity = expression(fifelse(
    is.na(p_HEIDI) & is.na(phet_ivw) & is.na(psigmay_mrlink2), NA,
    (p_HEIDI >= 0.05) %in% TRUE | (phet_ivw >= 0.05) %in% TRUE |
      (psigmay_mrlink2 >= 0.05) %in% TRUE
  )),
  coloc = expression(fifelse(
    is.na(coloc_pph4) & is.na(mvcoloc_pph4), NA,
    (coloc_pph4 >= 0.8) %in% TRUE | (mvcoloc_pph4 >= 0.8) %in% TRUE
  )),
  mr_sens = expression(mr & sensitivity),
  mr_coloc = expression(mr & coloc),
  mr_sens_coloc = expression(mr & sensitivity & coloc)
)

# Create columns for evidence criteria and count
for (e in names(evidence_criteria)) {
  df_msmr_tenk10k[, (e) := eval(evidence_criteria[[e]])]
}


# calculate max evidence
df_msmr_tenk10k[, `:=`(
    max_evidence = case_when(
      mr_sens_coloc ~ "mr_sens_coloc",
      mr_coloc      ~ "mr_coloc",
      mr_sens       ~ "mr_sens",
      mr            ~ "mr"
    ) |> factor(levels = c("mr", "mr_sens", "mr_coloc", "mr_sens_coloc"))
  )]

# Save intermediate results
write_parquet(df_msmr_tenk10k, snakemake@output[[1]], compression = "gzip")
