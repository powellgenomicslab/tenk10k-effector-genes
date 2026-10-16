# Trait-level insights from cell-type-specific effector gene mapping

These scripts look at the MR results across traits and genes. They compute trait–trait Spearman correlations of signed MR P-values, then count the effective number of independent traits per gene, which is a measure of pleiotropy. As an example, they examine the pleiotropic FADS/MYRF locus. The outputs are Figure 2 (panels a–b; panels c–d come from `scripts/3-pleiotropy/count_evidence.R`), Extended Data Figure 1, Supplementary Figure 7 and Supplementary Tables 4–5.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `fig_2_heatmap_spearman_calculation.R` | Trait–trait Spearman correlations (input to Fig 2a, Supp Tables 4–5) | intermediate | v5 parquet, `metadata/trait.tsv` | `results/rrho/spearman_corr_all_trait_combos_strictmr.tsv` |
| `fig_2_heatmap.R` | Fig 2a panel | intermediate | `preprocess_results.R` objects, Spearman table, `results/rrho/rrho_pvalues_all_trait_combos_two.sided.tsv`, `metadata/trait.tsv` | `figures/panels/2-trait/figure_2_heatmap_{disease_only,legend_disease}.rds` |
| `gene_features_pleiotropy.R` | Fig 2b panel; Extended Data Fig 1 panel | intermediate | v5 parquet, Spearman table, `metadata/` | `figures/panels/2-trait/{n_phenotypes_by_n_celltypes_split,fig_1_pleiotropic_genes_scatter}.rds` |
| `overview_combined.R` | Fig 2; Extended Data Fig 1 | `make fig-m2`, `make fig-e1` | panel RDS files above, plus `figures/panels/3-pleiotropy/evidence_tally_plot.rds` | `figures/fig-m2.{pdf,png}`, `figures/fig-e1.{pdf,png}` |
| `myrf_locus_supplementary_heatmap.R` | Supp Fig 7 | `make fig-s7` | v5 parquet, `metadata/` | `figures/fig-s7.{pdf,png}`, `figures/panels/2-trait/fads_locus_pleiotropic_scatter.{pdf,png}` |
| `tab_spearman_corr.R` | Supp Table 4 | `make tab-s4` | Spearman table; `metadata/table_column_names.xlsx` | `tables/` (Supp Table 4) |
| `tab_gene_pleiotropy.R` | Supp Table 5 | `make tab-s5` | `preprocess_results.R` objects, Spearman table; `metadata/table_column_names.xlsx` | `tables/` (Supp Table 5) |

## Run order

1. `fig_2_heatmap_spearman_calculation.R` writes the Spearman table that the other scripts read.
2. `fig_2_heatmap.R`, `gene_features_pleiotropy.R` and `scripts/3-pleiotropy/count_evidence.R` write the panel objects.
3. `overview_combined.R` assembles Figure 2 and Extended Data Figure 1.

`myrf_locus_supplementary_heatmap.R` is independent of the others. `tab_spearman_corr.R` and `tab_gene_pleiotropy.R` need only step 1.

## Restricted inputs

None. All scripts read only the aggregated summary-statistic results and `metadata/`.
