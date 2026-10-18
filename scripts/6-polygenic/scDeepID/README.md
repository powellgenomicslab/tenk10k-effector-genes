# Integration between scDRS and scDeepID

The scDeepID model is in a separate repository:
<https://github.com/powellgenomicslab/scDeepID_TenK10K_manuscripts>. Its B-cell cell
function scores are compared between SLE-enriched B cells (scDRS Monte Carlo
P < 0.05) and the other B cells sampled for scDRS. The test of each of the 400 cell
functions, with their GO top-level categories, comes from the scDeepID run as
`results/scDeepID/pathway_test_full.csv` and is in the data bundle.

| Script | Manuscript item | Make target | Key inputs | Output |
|---|---|---|---|---|
| `plot_scdeepid_categories.R` | Supp Fig 15 | `make fig-s15` | `results/scDeepID/pathway_test_full.csv` | `figures/fig-s15.{pdf,png}` |
| `plot_scdeepid_bcell.R` | Extended Data Fig 7 | `make fig-e7` (not yet: see below) | `results/scDeepID/tenk10k_bcell.sle.ed_fig7_pathways.anon.parquet`, `pathway_test_full.csv`, the scDRS cell-level files | `figures/fig-e7.{pdf,png}` |
| `tab_scdeep_id.R` | Supp Table 10 | `make tab-s10` | `results/scDeepID/pathway_test_full.csv`, `resources/enrichment/gmt/GO_bp_2024_terms.tsv` | `tables/tab-s10.tsv` |

Extended Data Fig 7 needs, for each sampled B cell, the scDeepID UMAP and the scores
of the five cell functions shown. Rule `scdeepid_bcell_extract`
(`workflow/rules/scdrs.smk`) writes them from the scDeepID scores, which are not
shared, and checks its P values against `pathway_test_full.csv` to make sure the
scores come from the same run. Until that file is in the bundle, `fig_e7` is listed
in `NOT_REBUILDABLE` in `scripts/Snakefile`.

The figures were first drawn with matplotlib by
`scripts/4-polygenic/scDeepID/scDRS_scDeepID_BCELL.py` (`main` branch).
