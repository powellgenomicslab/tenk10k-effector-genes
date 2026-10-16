# Polygenic enrichment reveals trait-relevant immune cell types and gene programs

These scripts turn the single-cell Disease Relevance Score (scDRS) output of the
pipeline (`workflow/rules/scdrs.smk`) into the manuscript's polygenic enrichment
figures and tables. They cover enrichment of the 100 traits across 28 immune cell
types in TenK10K, the Crohn's disease example, two robustness checks (an
alternative cell-sampling strategy and replication in the AIFI Immune Health
Atlas), and the supplementary table of scDeepID cell functions in SLE-enriched
B cells.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `scDRS/plot_scdrs_overall.R` | Fig 4a–c; Extended Data Fig 6; Supp Fig 14; Supp Table 9 | `make fig-m4`, `make fig-e6`, `make fig-s14`, `make tab-s9` | `results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv`, `results/aggregate/tenk10k_phase1.scdrs.cell_score.anon.parquet.gz`, `results/aggregate/tenk10k_phase1.scdrs.cell_mcp.anon.parquet.gz`, `metadata/{cell,trait,trait_category}.tsv` | `figures/fig-m4.{pdf,png}`, `figures/fig-e6.{pdf,png}`, `figures/fig-s14.{pdf,png}`, `tables/tab-s9.tsv` |
| `scDRS/compare_scdrs_aifi.R` | Supp Fig 23 (cells sampled per cell type); Supp Fig 24 (scDRS Z-score, TenK10K vs AIFI) | `make fig-s23`, `make fig-s24` | `results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv`, `results/aggregate/immune_health_atlas_annotated.scdrs.cell_type_stats.tsv`, `metadata/{cell,trait,trait_category}.tsv` | `figures/fig-s23.{pdf,png}`, `figures/fig-s24.{pdf,png}` |
| `scDRS/aifi-labelTransfer.py` | Input to Supp Figs 23–24 (CellTypist label transfer onto AIFI) | (intermediate) | `data/immune_health_atlas.h5ad`, `output/wg2_scpred_model.pkl` | `output/immune_health_atlas_annotated.h5ad` |
| `scDRS/scdrs_supp_sampling.R` | Supp Figs 21–22 (Crohn's disease scDRS across two sampling strategies) | `make fig-s21`, `make fig-s22` | per-cell scDRS scores and MC P-values, `results/misc/tenk10k_prop_sample.crohns.cell_score.anon.parquet.gz`, `metadata/cell.tsv` | `figures/fig-s21.{pdf,png}`, `figures/fig-s22.{pdf,png}` |
| `scDeepID/tab_scdeep_id.R` | Supp Table 10 (scDeepID pathway scores in SLE-enriched B cells) | `make tab-s10` | `results/scDeepID/pathway_test_full.csv`, `resources/enrichment/gmt/GO_bp_2024_terms.tsv` | `tables/tab-s10.tsv` |

Extended Data Fig 7a–b and Supp Fig 15 (scDeepID cell function scores) are
produced in the scDeepID repository, not here; see `scDeepID/README.md`.

All `write_table()` calls rename columns from `metadata/table_column_names.xlsx`.

## Run order

1. The pipeline computes scDRS for TenK10K (`tenk10k_phase1`, and
   `tenk10k_prop_sample` for the sampling check).
2. `scDRS/plot_scdrs_overall.R` and `scDRS/scdrs_supp_sampling.R`.
3. For the AIFI replication: run `scDRS/aifi-labelTransfer.py`, then the
   pipeline's scDRS on `immune_health_atlas_annotated`, then
   `scDRS/compare_scdrs_aifi.R`.
4. `scDeepID/tab_scdeep_id.R` needs the scDeepID output
   (`results/scDeepID/pathway_test_full.csv`) from the external repository.

## Restricted inputs

- `plot_scdrs_overall.R` and `scdrs_supp_sampling.R` read **per-cell** scDRS
  scores and Monte Carlo P-values. These are keyed to individual TenK10K donors
  and are not part of the public release, so a third party cannot run these two
  scripts. Only the cell-type-level statistics are released.
- `aifi-labelTransfer.py` reads the AIFI atlas `.h5ad` (from its own source) and a
  CellTypist model trained on TenK10K single-cell data.
