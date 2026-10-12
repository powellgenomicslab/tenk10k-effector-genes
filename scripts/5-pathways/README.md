# Identifying disease-relevant biological pathways using single-cell genetics

These scripts summarise the rank-based STRING functional enrichment of MR results in CD4+ cytotoxic T lymphocytes (CD4 CTLs) for 69 diseases. A bipartite network links the non-HLA enriched terms to diseases. The STRING protein–protein interaction network, with housekeeping genes removed, highlights the genes in multiple sclerosis-enriched terms and compares the betweenness centrality of MR-nominated genes with that of other genes (Welch's two-sample t-test).

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `enrichment_features.R` | Extended Data Fig. 5; Supplementary Figs 12 and 13 | `make fig-e5`, `make fig-s12`, `make fig-s13` | `results/preprocessed/tenk10k_phase1.v5.parquet.gz`; `metadata/{trait.tsv,cell.tsv,trait_category.tsv,gencode.v44.gene_type.tsv.gz}`; `results/rrho/spearman_corr_all_trait_combos_strictmr.tsv`; `results/enrichment/tenk10k_phase1_mr_strict/mr_strict/CD4_CTL.string_rank_api{.tsv,.ppi_network.tsv,.bipartite.tsv}`; `resources/enrichment/string_id_map.tsv`; `resources/misc/HOUNKPE_HOUSEKEEPING_GENES.v2026.1.Hs.tsv` | `figures/fig-e5.{png,pdf}`, `figures/fig-s12.{png,pdf}`, `figures/fig-s13.{png,pdf}` |
| `tab_enrich_cd4ctl.R` | Supplementary Table 8 | `make tab-s8` | `metadata/trait.tsv`; `results/enrichment/tenk10k_phase1_mr_strict/mr_strict/CD4_CTL.string_rank_api.tsv`; `metadata/table_column_names.xlsx` | `tables/tab-s8.tsv` |

## Run order

`enrichment_features.R` reads `results/rrho/spearman_corr_all_trait_combos_strictmr.tsv`, which `scripts/2-trait/fig_2_heatmap_spearman_calculation.R` writes, so run that script first. `tab_enrich_cd4ctl.R` is independent.

## Restricted inputs

None. Both scripts read summary-level results only.
