# Polygenic enrichment reveals trait-relevant immune cell types and gene programs

These scripts turn the single-cell Disease Relevance Score (scDRS) output of the
pipeline (`workflow/rules/scdrs.smk`) into the manuscript's polygenic enrichment
figures and tables. They cover enrichment of the 100 traits across 28 immune cell
types in TenK10K, the Crohn's disease example, two robustness checks (an
alternative cell-sampling strategy and replication in the AIFI Immune Health
Atlas), and scDeepID cell functions in SLE-enriched B cells.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `scDRS/plot_scdrs_overall.R` | Fig 4a–c; Extended Data Fig 6; Supp Fig 14; Supp Table 9 | `make fig-m4`, `make fig-e6`, `make fig-s14`, `make tab-s9` | `results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv`, `results/aggregate/tenk10k_phase1.scdrs.cell_score.anon.parquet.gz`, `results/aggregate/tenk10k_phase1.scdrs.cell_mcp.anon.parquet.gz`, `metadata/{cell,trait,trait_category}.tsv` | `figures/fig-m4.{pdf,png}`, `figures/fig-e6.{pdf,png}`, `figures/fig-s14.{pdf,png}`, `tables/tab-s9.tsv` |
| `scDRS/compare_scdrs_aifi.R` | Supp Fig 23 (cells sampled per cell type); Supp Fig 24 (scDRS Z-score, TenK10K vs AIFI) | `make fig-s23`, `make fig-s24` | `results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv`, `results/aggregate/immune_health_atlas_annotated.scdrs.cell_type_stats.tsv`, `metadata/{cell,trait,trait_category}.tsv` | `figures/fig-s23.{pdf,png}`, `figures/fig-s24.{pdf,png}` |
| `workflow/rules/snakescripts/scdrs/aifi-labelTransfer.py` (analysis pipeline) | Input to Supp Figs 23–24 (CellTypist label transfer onto AIFI) | (intermediate) | `data/immune_health_atlas.h5ad`, `output/wg2_scpred_model.pkl` | `output/immune_health_atlas_annotated.h5ad` |
| `scDRS/scdrs_supp_sampling.R` | Supp Figs 21–22 (Crohn's disease scDRS across two sampling strategies) | `make fig-s21`, `make fig-s22` | per-cell scDRS scores and MC P-values, `results/misc/tenk10k_prop_sample.crohns.cell_score.anon.parquet.gz`, `metadata/cell.tsv` | `figures/fig-s21.{pdf,png}`, `figures/fig-s22.{pdf,png}` |
| `scDeepID/plot_scdeepid_categories.R` | Supp Fig 15 (scDeepID cell function scores by GO category in SLE-enriched B cells) | `make fig-s15` | `results/scDeepID/pathway_test_full.csv` | `figures/fig-s15.{pdf,png}` |
| `scDeepID/plot_scdeepid_bcell.R` | Extended Data Fig 7 (scDeepID cell function scores in SLE-enriched B cells) | `make fig-e7` | `results/scDeepID/{tenk10k_bcell.sle.ed_fig7_pathways.anon.parquet,pathway_test_full.csv}`, the scDRS cell-level files | `figures/fig-e7.{pdf,png}` |
| `scDeepID/tab_scdeep_id.R` | Supp Table 10 (scDeepID pathway scores in SLE-enriched B cells) | `make tab-s10` | `results/scDeepID/pathway_test_full.csv`, `resources/enrichment/gmt/GO_bp_2024_terms.tsv` | `tables/tab-s10.tsv` |

The scDeepID inputs are described in `scDeepID/README.md`.

All `write_table()` calls rename columns from `metadata/table_column_names.xlsx`.

## Run order

1. The pipeline computes scDRS for TenK10K (`tenk10k_phase1`, and
   `tenk10k_prop_sample` for the sampling check).
2. `scDRS/plot_scdrs_overall.R` and `scDRS/scdrs_supp_sampling.R`.
3. For the AIFI replication: run `workflow/rules/snakescripts/scdrs/aifi-labelTransfer.py`, then the
   pipeline's scDRS on `immune_health_atlas_annotated`, then
   `scDRS/compare_scdrs_aifi.R`.
4. For scDeepID: `scDeepID/tab_scdeep_id.R` and `scDeepID/plot_scdeepid_categories.R`
   read `results/scDeepID/pathway_test_full.csv`, written by the scDeepID revision notebook;
   `scDeepID/plot_scdeepid_bcell.R` reads the output of the pipeline's
   `scdeepid_bcell_extract`.

## Restricted inputs

- `plot_scdrs_overall.R` and `scdrs_supp_sampling.R` read **per-cell** scDRS
  scores, Monte Carlo P-values and UMAP coordinates. The pipeline computes them
  (`workflow/rules/scdrs.smk`), and they are shared in the data bundle with each
  cell ID made of the 10x barcode and an anonymous code in place of the
  sequencing pool (rule `scdrs_anonymise_cell_ids`).
- `aifi-labelTransfer.py` (in the pipeline, `workflow/rules/snakescripts/scdrs/`)
  reads the AIFI atlas `.h5ad` (from its own source) and a CellTypist model
  trained on TenK10K single-cell data; neither is shared.
