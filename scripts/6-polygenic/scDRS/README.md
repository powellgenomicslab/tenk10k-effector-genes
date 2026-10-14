# scDRS: polygenic enrichment at single-cell resolution

These scripts build the manuscript's scDRS figures and tables from the
pipeline output. The full table of manuscript items is in `../README.md`.

| Script | Purpose |
|---|---|
| `plot_scdrs_overall.R` | Fig 4, Extended Data Fig 6, Supp Fig 14 and Supp Table 9 |
| `compare_scdrs_aifi.R` | Replication in the AIFI Immune Health Atlas (Supp Figs 23–24) |
| `aifi-labelTransfer.py` | Transfers TenK10K cell-type labels onto the AIFI atlas. It writes the annotated `.h5ad` that the pipeline scores and `compare_scdrs_aifi.R` then summarises |
| `scdrs_supp_sampling.R` | Crohn's disease robustness to the cell-sampling strategy |

The pipeline computes the scores themselves, not these scripts. See
`workflow/rules/scdrs.smk` and `workflow/rules/snakescripts/scdrs/`. Per-dataset
parameters are in `resources/scdrs/config/{study}.yaml`.

> [!IMPORTANT]
> `plot_scdrs_overall.R` and `scdrs_supp_sampling.R` read **per-cell** scDRS
> scores. These are keyed to individual TenK10K donors and are **not** part of
> the public data release, so a third party cannot run these two scripts. Only
> the cell-type-level statistics are released. `aifi-labelTransfer.py` needs the
> AIFI atlas `.h5ad`, which you obtain from its own source.
