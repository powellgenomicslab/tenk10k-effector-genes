# supplementary table: per-gene pleiotropy summary (n cell types / phenotypes / effective
# phenotypes each MR gene is associated with)
#
# Minimal gene_summary build shared with scripts/2-trait/gene_features_pleiotropy.R
# and scripts/5-pathways/enrichment_features.R.
# Sourced via preprocess_results.R (not a second hardcoded parquet path) so this stays on
# the same MR result set (currently v5) as tab_mr_overlap.R / tab_smr_overview.R.

source("scripts/0-preprocess/preprocess_results.R")
source("scripts/util/write_table.R")

# phenotype correlation matrix, used to count effective (non-redundant) phenotypes per gene
df_rrho <- read_tsv("results/rrho/spearman_corr_all_trait_combos_strictmr.tsv") %>%
  select(trait1 = p1, trait2 = p2, value = spearman_corr)

rrho_traits <- union(df_rrho$trait1, df_rrho$trait2)
mat_cor <- bind_rows(
    df_rrho,
    rename(df_rrho, trait1 = trait2, trait2 = trait1),
    tibble(trait1 = rrho_traits, trait2 = rrho_traits, value = 1)
  ) %>%
  pivot_wider(names_from = trait2, values_from = value) %>%
  column_to_rownames("trait1") %>%
  data.matrix() %>%
  .[rrho_traits, rrho_traits]

# Li and Ji (2005) effective number of independent tests, via eigenvalue decomposition
count_effective_tests <- function(pheno_vec, cor_matrix) {
  if (length(pheno_vec) <= 1) return(length(pheno_vec))
  pheno_in_matrix <- intersect(pheno_vec, rownames(cor_matrix))
  if (length(pheno_in_matrix) <= 1) return(length(pheno_in_matrix))
  sub_cor <- cor_matrix[pheno_in_matrix, pheno_in_matrix]
  eigenvalues <- eigen(sub_cor, symmetric = TRUE, only.values = TRUE)$values
  f <- function(x) {
    if (x <= 0) return(0)
    indicator <- as.numeric(x >= 1)
    fractional <- x - floor(x)
    indicator + fractional
  }
  sum(sapply(eigenvalues, f))
}

gene_summary <- df_msmr_tenk10k[mr == TRUE, .(
  n_phenotypes = length(unique(pheno_label)),
  n_pheno_cat = length(unique(pheno_cat)),
  n_celltypes = length(unique(biosample)),
  phenotypes = list(unique(pheno_label)),
  celltypes = list(unique(biosample))
), by = probeID]

gene_summary[, n_effective_phenotypes := sapply(phenotypes, function(p) {
  count_effective_tests(p, mat_cor)
})]

gene_summary[df_gene_annot, gene := i.hgnc_symbol, on = c("probeID" = "ensembl_gene_id")]

tab_gene_pleiotropy <- gene_summary %>%
  mutate(
    phenotypes = sapply(phenotypes, paste, collapse = "; "),
    celltypes = sapply(celltypes, paste, collapse = "; ")
  ) %>%
  as.data.table()

write_table(tab_gene_pleiotropy, "gene_pleiotropy", 5)
