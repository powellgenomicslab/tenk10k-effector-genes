# scDRS: polygenic enrichment at single-cell resolution

These scripts build the manuscript's scDRS figures and tables from the
pipeline output. The full table of manuscript items is in `../README.md`.

| Script | Purpose |
|---|---|
| `plot_scdrs_overall.R` | Fig 4, Extended Data Fig 6, Supp Fig 14 and Supp Table 9 |
| `compare_scdrs_aifi.R` | Replication in the AIFI Immune Health Atlas (Supp Figs 23–24) |
| `workflow/rules/snakescripts/scdrs/aifi-labelTransfer.py` (analysis pipeline) | Transfers TenK10K cell-type labels onto the AIFI atlas. It writes the annotated `.h5ad` that the pipeline scores and `compare_scdrs_aifi.R` then summarises |
| `scdrs_supp_sampling.R` | Crohn's disease robustness to the cell-sampling strategy |

The pipeline computes the scores themselves, not these scripts. See
`workflow/rules/scdrs.smk` and `workflow/rules/snakescripts/scdrs/`. Per-dataset
parameters are in `resources/scdrs/config/{study}.yaml`.

> [!IMPORTANT]
> `plot_scdrs_overall.R` and `scdrs_supp_sampling.R` read **per-cell** scDRS
> scores and UMAP coordinates computed by the pipeline. They are shared in the
> data bundle with anonymised cell IDs (10x barcode and an anonymous code in place
> of the sequencing pool). The label transfer onto AIFI is a pipeline step
> (`workflow/rules/snakescripts/scdrs/aifi-labelTransfer.py`) and needs the AIFI
> atlas `.h5ad`, which you obtain from its own source.
