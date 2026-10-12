# supplementary tables: TenK10K MR vs. MAGMA (GWAS) and vs. eQTLgen bulk-blood MR
#
# Recovered from scripts/mr_overview/overview_supplement.R (deleted in commit e137f97,
# 2026-05-12), which wrote this analysis to a single combined sheet "mr_gwas_eqtlgen".
# The current column-label template (metadata/table_column_names.xlsx) splits
# it into three sheets instead: mr_gwas_overlap, mr_magma, mr_eqtlgen.
#
# mr_magma (tbl_order 6) and mr_eqtlgen (tbl_order 7) were confirmed against the live
# spreadsheet's existing "6-mr_magma"/"7-mr_eqtlgen" tabs: identical trait set and gene
# counts for all 99 traits present in both. (The live tabs show one trait as "PCOS";
# this table correctly shows "PMOS" per the current trait metadata label --
# the live sheet is just stale on that one label and will pick up "PMOS" when re-written.)
#
# mr_gwas_overlap is not a manuscript Supplementary Table and is not written here.

source("scripts/0-preprocess/preprocess_results.R")
source("scripts/util/write_table.R")

# --- mr_magma / mr_eqtlgen: per-trait gene tally split by MAGMA / eQTLgen concordance ---

df_gene_summary_by_pheno <- df_msmr_tenk10k[
  probeID %in% gene_universe,
  .(mr = max(mr),
    magma_gene = max(magma_gene),
    eqtlgen_mr = max(eqtlgen_mr)),
  by = list(probeID, phenotype, pheno_label, pheno_cat)
] |>
  mutate(
    eqtlgen = case_when(
      eqtlgen_mr == 1 & mr == 1 ~ "TenK10K MR & eQTLgen MR",
      eqtlgen_mr == 0 & mr == 1 ~ "TenK10K MR only",
      eqtlgen_mr == 1 & mr == 0 ~ "eQTLgen MR only",
      TRUE ~ NA_character_),
    gwas = case_when(
      magma_gene == 1 & mr == 1 ~ "TenK10K MR & MAGMA",
      magma_gene == 0 & mr == 1 ~ "TenK10K MR only",
      magma_gene == 1 & mr == 0 ~ "MAGMA only",
      TRUE ~ NA_character_)
  ) |>
  pivot_longer(c(gwas, eqtlgen),
               values_drop_na = TRUE,
               names_to = "annotation", values_to = "value") |>
  group_by(phenotype, pheno_label, pheno_cat, annotation, value) |>
  tally(name = "n_genes") |>
  ungroup() |>
  complete(nesting(phenotype, pheno_label, pheno_cat),
           nesting(annotation, value),
           fill = list(n_genes = 0)) |>
  group_by(phenotype, pheno_label, pheno_cat, annotation) |>
  mutate(n_genes_total = sum(n_genes),
         prop_genes = n_genes / sum(n_genes)) |>
  ungroup() |>
  mutate(fct_pheno = fct_reorder(pheno_label, n_genes_total, mean) |> fct_rev()) |>
  setDT()

# id_cols must be pinned explicitly: n_genes_total differs between the gwas and eqtlgen
# sub-groups of the same trait, so leaving pivot_wider() to infer id columns from
# "everything else" produces two half-NA rows per trait instead of one full row.
tab_gene_summary_by_pheno <- df_gene_summary_by_pheno %>%
  pivot_wider(id_cols = c(fct_pheno, pheno_cat, pheno_label, phenotype),
              names_from = c(annotation, value),
              values_from = c(n_genes, prop_genes),
              names_sep = ".") %>%
  setDT()

write_table(tab_gene_summary_by_pheno, "mr_magma", 6)
write_table(tab_gene_summary_by_pheno, "mr_eqtlgen", 7)
