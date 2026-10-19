# Integration between scDRS and scDeepID

The scDeepID model is in a separate repository:
<https://github.com/powellgenomicslab/scDeepID_TenK10K_manuscripts>. Its B-cell cell
function scores (run `exp_10k_tr_revision/bcell_human_gobp_2024_ep10_lr0.02_bs32_rentention0.0_maxgs400`)
are compared between SLE-enriched B cells (scDRS Monte Carlo P < 0.05) and the other
B cells sampled for scDRS, by the scDeepID revision notebook
(`revision_cleanup_v3_github.ipynb` in the run directory). The notebook writes
`results/scDeepID/pathway_test_full.csv`: a one-sided t-test, BH FDR and group means
for each of the 400 cell functions, and their GO categories (the children of
"biological process", GO:0008150, with at least 5 cell functions at FDR < 0.05, with
membership taken through the GO graph in `go-basic.obo`). That file is in the data
bundle.

| Script | Manuscript item | Make target | Key inputs | Output |
|---|---|---|---|---|
| `plot_scdeepid_categories.R` | Supp Fig 15 | `make fig-s15` | `results/scDeepID/pathway_test_full.csv` | `figures/fig-s15.{pdf,png}` |
| `plot_scdeepid_bcell.R` | Extended Data Fig 7 | `make fig-e7` | `results/scDeepID/tenk10k_bcell.sle.ed_fig7_pathways.anon.parquet`, `pathway_test_full.csv`, the scDRS cell-level files | `figures/fig-e7.{pdf,png}` |
| `tab_scdeep_id.R` | Supp Table 10 | `make tab-s10` | `results/scDeepID/pathway_test_full.csv`, `resources/enrichment/gmt/GO_bp_2024_terms.tsv` | `tables/tab-s10.tsv` |

Supp Fig 15 shows, for each category, the 30 cell functions with the lowest FDR
(FDR < 0.05). Extended Data Fig 7 shows the five with the lowest FDR in "biological
regulation". It needs, for each sampled B cell, the scDeepID UMAP and the scores of
those five; rule `scdeepid_bcell_extract` (`workflow/rules/scdrs.smk`) writes them
from the scDeepID scores, which are not shared, and checks its P values against
`pathway_test_full.csv`.

The figures were drawn with matplotlib by the notebook; these scripts redraw them
with ggplot2 from the same values.
