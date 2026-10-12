# scDRS: polygenic enrichment at single-cell resolution

Scripts that produce the manuscript's scDRS figures and tables from pipeline
output.

| Script | Purpose |
|---|---|
| `scdrs_main_supp.R` | Main and supplementary scDRS figures |
| `scdrs_supp_sampling.R` | Sensitivity of the cell-type results to cell sampling |
| `compare_scdrs_aifi.R` | Replication against the AIFI Immune Health Atlas |
| `aifi-labelTransfer.py` | Transfers TenK10K cell-type labels onto the AIFI atlas, producing the annotated `.h5ad` that `compare_scdrs_aifi.R` consumes |

The scores themselves are computed by the pipeline, not by these scripts: see
`workflow/rules/scdrs.smk` and `workflow/rules/snakescripts/scdrs/`, with
per-dataset parameters in `resources/scdrs/config/{study}.yaml`.

> [!IMPORTANT]
> `scdrs_main_supp.R` and `scdrs_supp_sampling.R` read **per-cell** scores
> (`results/aggregate/tenk10k_phase1.scdrs.cell_score.tsv.parquet.gz` and
> `.cell_mcp.tsv.parquet.gz`). Those are keyed to individual TenK10K donors and
> are **not** part of the public data release, so these two scripts cannot be run
> by a third party. Only the aggregated cell-type-level statistics are released.
> `aifi-labelTransfer.py` likewise needs the AIFI atlas `.h5ad`, obtained from
> its own source.
