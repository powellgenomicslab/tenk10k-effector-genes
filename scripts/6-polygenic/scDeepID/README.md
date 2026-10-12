# Integration between scDRS and scDeepID

The scDeepID model and the cell function scores for SLE-enriched B cells are
produced in a separate repository:
<https://github.com/powellgenomicslab/scDeepID_TenK10K_manuscripts>. That
repository also makes the figures: Extended Data Fig 7a–b and Supp Fig 15.

| Script | Manuscript item | Make target | Key inputs | Output |
|---|---|---|---|---|
| `tab_scdeep_id.R` | Supp Table 10 | `make tab-s10` | `results/scDeepID/pathway_test_full.csv` (from scDeepID), `resources/enrichment/gmt/GO_bp_2024_terms.tsv` | `tables/tab-s10.tsv` |
