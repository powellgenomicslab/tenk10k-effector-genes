# Assessing horizontal pleiotropy in single-cell eQTL Mendelian randomisation via sensitivity analyses and colocalisation

These scripts summarise how far the MR associations are supported by the three sensitivity analyses (HEIDI, Cochran's Q and MR-link-2) and by colocalisation (coloc and coloc-SuSiE). They also show how the overlap between MR and colocalisation changes with the significance thresholds, and how the number of instruments affects MR and sensitivity statistics. The evidence flags (`mr`, `mr_sens`, `mr_sens_coloc`, `max_evidence` and so on) are read directly from the preprocessed v5 parquet.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `count_evidence.R` | Fig. 2c–d (evidence tally panel); Extended Data Fig. 2 (support by trait category) | `make fig-e2` (the panel is consumed by `make fig-m2`) | `results/preprocessed/tenk10k_phase1.v5.parquet.gz` | `figures/panels/3-pleiotropy/evidence_tally_plot.rds`, `figures/fig-e2.{png,pdf}` |
| `coloc_mr_overlap.R` | Supplementary Fig. 8a–c | `make fig-s8` | `scripts/0-preprocess/preprocess_results.R` inputs | `figures/fig-s8.{png,pdf}` |
| `mr_stat_by_n_instrument.R` | Extended Data Fig. 3a–d | `make fig-e3` | `scripts/0-preprocess/preprocess_results.R` inputs | `figures/fig-e3.{png,pdf}` |

## Run order

`count_evidence.R` must run before `scripts/2-trait/overview_combined.R`, which reads `figures/panels/3-pleiotropy/evidence_tally_plot.rds` as Fig. 2c. The other two scripts are independent.

## Restricted inputs

None. All three scripts read summary-level results only.
