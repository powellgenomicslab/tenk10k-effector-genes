# preprocessing script with stringent MSMR results as the baseline results

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(readxl)
library(qvalue)

setDTthreads(Sys.getenv("NCPUS"))
df_cell_map <- fread("resources/metadata/cell_map.tsv")
df_trait_map_all <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

setDT(df_trait_map)
phenotypes <- df_trait_map$trait_id

cat_order <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx",
  sheet = "trait_category_order"
) %>%
  pull(cat_order)

df_msmr_tenk10k <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |> 
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

# mrlink2
df_mrlink2 <- read_parquet("results/aggregate/tenk10k_phase1_sensitivity.mrlink2.parquet.gz")

# ivw-ld
df_ivw <- read_parquet("/g/data/fy54/analysis/tenk10k-causal/results/aggregate/tenk10k_phase1_sensitivity.ivw-ld.parquet.gz")

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
  setDT()

# annotate ivw-ld results
df_msmr_tenk10k[df_ivw, `:=`(b_ivw = i.estimate, se_ivw = i.se, p_ivw = i.pval, phet_ivw = i.het_pval),
  on = c("biosample", "phenotype", "probeID")]

# annotate mrlink2 results
df_msmr_tenk10k[df_mrlink2, `:=`(
    b_mrlink2 = i.alpha,
    se_mrlink2 = `i.se(alpha)`,
    p_mrlink2 = `i.p(alpha)`,
    psigmay_mrlink2 = `i.p(sigma_y)`
  ), on = c("biosample", "phenotype", "probeID")]

# annotate coloc  & mv coloc results
df_msmr_tenk10k[df_coloc, coloc_pph4 := i.PP.H4.abf, on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]
df_msmr_tenk10k[df_mvcoloc, mvcoloc_pph4 := i.PP.H4.abf, on = c("probeID" = "gene", "phenotype" = "pheno", "biosample")]

# arrange and rename columns

cols <- c("biosample", "cell_type", "major_cell_type",
          "phenotype", "pheno_label", "pheno_cat", "supercategory",
          "probeID", "Gene", "gene_type", "ProbeChr", "Probe_bp",
          "topSNP", "topSNP_chr", "topSNP_bp",
          "A1", "A2", "Freq",
          "b_GWAS", "se_GWAS", "p_GWAS",
          "b_eQTL", "se_eQTL", "p_eQTL",
          "b_SMR", "se_SMR", "p_SMR", "p_SMR_multi", "p_HEIDI", "nsnp_HEIDI",
          "qval_msmr_pheno", "lfdr_msmr_pheno",
          "b_ivw", "se_ivw", "p_ivw", "phet_ivw",
          "b_mrlink2", "se_mrlink2", "p_mrlink2", "psigmay_mrlink2",
          "coloc_pph4", "mvcoloc_pph4")

df_msmr_tenk10k <- df_msmr_tenk10k[, ..cols]

evidence_criteria <- list(
  mr = expression(lfdr_msmr_pheno < 0.05),
  sensitivity = expression(p_HEIDI >= 0.05 | phet_ivw >= 0.05 | psigmay_mrlink2 >= 0.05),
  coloc = expression(coloc_pph4 >= 0.8 | mvcoloc_pph4 >= 0.8),
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
