# TenK10K Phase1 causal inference manuscript

![Static Badge](https://img.shields.io/badge/version-preprint_v4-red)
[![medRxiv](https://img.shields.io/badge/medRxiv-10.1101%2F2025.08.28.25334614-B31B1B)](https://www.medrxiv.org/content/10.1101/2025.08.28.25334614v4)

> [!NOTE]
> This repository is still under active development and will be updated from time to time

## Study design
![](figures/biorender/study_design.png)

## Overview

This repository contains code, data, and workflows for the TenK10K causal inference manuscript. The project is organized into several directories, each corresponding to a major analysis step or component.

## Directory Structure

- **`figures/`**  
  Contains figures generated for the manuscript, including study design diagrams.

- **`metadata/`**  
  Contains metadata files for traits, trait categories, and cell types from TenK10K Phase 1 sc-eQTL analysis.

- **`scripts/`**  
  Main analysis scripts, organized by analysis section:
  - **`0-preprocess/`** - Shared preprocessing of the aggregate results
  - **`1-overview/`** - Study overview and summary statistics  
  - **`2-mr/`** - Mendelian randomisation analyses and comparisons
  - **`3-comparison/`** - Overlap between MR, coloc, MAGMA and eQTLGen
  - **`4-polygenic/`** - Polygenic enrichment analyses (scDRS, scDeepID)
  - **`5-drug/`** - Drug target enrichment and therapeutic relevance
  - **`6-crohns/`** - Crohn's disease case study with matched single-cell data
  - **`util/`** - Utility functions and helper scripts

- **`metadata/`**, **`config/`**, **`workflow/envs/`**, **`resources/{misc,scdrs/config}/`**  
  Small tracked configuration and metadata that the code needs in order to run:
  trait and cell-type definitions, gene annotation, conda environment specs,
  file-map templates and analysis parameters.

- **`sensitivity`**
  Sensitivity Snakemake pipeline for running [IVW-MR](https://mrcieu.github.io/TwoSampleMR/index.html) and [MR-link-2](https://github.com/adriaan-vd-graaf/mrlink2). This workflow reads in the intermediate results of the main MR pipeline to extract genes and instruments sets. sc-eQTL and GWAS summary statistics are also reformatted for use with sensitivity MR methods. 
    
- **`workflow/`**  
  Snakemake pipeline for reproducible data processing and analysis. The workflow handles data formatting, quality control, statistical analyses, and intermediate file generation.

## Usage

### Working directory

Every path in this repository is relative to a **working directory** that also
holds `resources/` (inputs) and `results/` (outputs). A bare checkout is not that
directory. Create one, place the released summary statistics and your own input
data under `resources/` and `results/`, and run everything from there. Both
pipelines are invoked from the working-directory root, not from inside
`workflow/`.

### Quick Start

1. **Run the Snakemake workflow** (generates the aggregate results). Requires
   Snakemake ≥ 8; see [`workflow/README.md`](workflow/README.md) for the
   prerequisites you must supply first.
   ```bash
   snakemake --snakefile workflow/snakefile \
             --profile workflow/profiles/default \
             --cores 8
   ```

2. **Run analysis scripts** (generates figures and tables):
   ```bash
   # Shared preprocessing - sourced by nearly every script below
   Rscript scripts/0-preprocess/preprocess_results.R

   # Overview analysis
   Rscript scripts/1-overview/study_design.R

   # Mendelian randomisation
   Rscript scripts/2-mr/mr_results_main.R

   # Method comparison
   Rscript scripts/3-comparison/coloc_mr_overlap.R

   # Polygenic analyses
   Rscript scripts/4-polygenic/scDRS/scdrs_main_supp.R

   # Drug target analysis
   Rscript scripts/5-drug/otp_combined.R

   # Crohn's disease case study
   Rscript scripts/6-crohns/crohns_case_study/plot_figures/3-annotated_heatmap.R
   ```

### Two version numbers

These are independent and easily confused:

- **Code / manuscript version** — tracked by this repository's branches. The
  current manuscript code is `v3`.
- **Summary-statistics version** — the vintage of the integrated results table,
  currently `v5` (`results/preprocessed/tenk10k_phase1.v5.parquet.gz`). Its
  provenance and the changes from `v4` are documented in
  [`CHANGELOG.v5.md`](CHANGELOG.v5.md). All analysis scripts read `v5`.



### Detailed Instructions

For detailed instructions on each analysis step, see the README files in each subdirectory:
- [`workflow/README.md`](workflow/README.md) - Snakemake pipeline setup and execution
- [`scripts/0-preprocess/README.md`](scripts/0-preprocess/README.md) - Data preprocessing
- [`scripts/1-overview/README.md`](scripts/1-overview/README.md) - Study overview and design
- [`scripts/2-mr/README.md`](scripts/2-mr/README.md) - Mendelian Randomization analysis
- [`scripts/4-polygenic/README.md`](scripts/4-polygenic/README.md) - Polygenic enrichment analysis
- [`scripts/5-drug/README.md`](scripts/5-drug/README.md) - Drug target enrichment
- [`scripts/6-crohns/README.md`](scripts/6-crohns/README.md) - Crohn's disease case study
- [`sensitivity/README.md`](sensitivity/README.md) - Mendelian Randomization - Sensitivity analysis


## Data Availability

Aggregated summary statistics — the MR/SMR, colocalisation and related result
tables — will be deposited on Zenodo and/or Hugging Face following manuscript
publication. The sc-eQTL summary statistics are released separately by the
mapping study (Cuomo et al.).

### What is not included, and why

This repository is code only, and the data release is deliberately limited to
aggregate statistics. The following are therefore **not** available, and the
scripts that depend on them cannot be run by a third party:

| Not released | Consequence |
|---|---|
| Individual-level single-cell expression (`.h5ad`) and TenK10K genotypes | The scDRS score computation and the locus-zoom scripts cannot be re-run |
| **Per-cell** scDRS scores (`*.scdrs.cell_score`, `*.scdrs.cell_mcp`) — keyed to individual donors | `scripts/4-polygenic/scDRS/scdrs_main_supp.R` and `scdrs_supp_sampling.R` cannot be run; only the aggregated cell-type statistics are released |
| Donor-derived LD matrices and per-locus dumps | `Crohns_example_locus_zoom.R`, `5-GPX1_locus_zoom.R` and the GWAS locus panels cannot be re-run |
| Third-party single-cell atlases (Kong et al. 2023 colon, AIFI Immune Health Atlas) | Obtain from their original sources |

A small number of curated files are still missing from this repository and must
be obtained separately or supplied by the user:

- `metadata/table_column_names.xlsx` — column-label mapping used by
  `scripts/util/write_table.R`; without it, supplementary-table generation fails.
- `metadata/otp.xlsx` — therapeutic-area and clinical-phase mappings for the
  drug analysis.
- `resources/misc/chembl_action_type.tsv`, and
  `results/otp/26.03/otp_concordance_by_phase.tsv` (produced by a commented-out
  block in `scripts/5-drug/otp_combined.R`).
- Public downloads to fetch yourself: the GENCODE v44 GTF, MSigDB
  `HOUNKPE_HOUSEKEEPING_GENES`, gnomAD v4.1 constraint metrics, and the
  UCSC hg19→hg38 liftover chain.

### Known gaps in the pipeline graph

Stated plainly, so that nobody assumes otherwise:

- **`resources/ma/` is a curated input**, not a pipeline output. See
  [`scripts/0-preprocess/gwas/README.md`](scripts/0-preprocess/gwas/README.md)
  for the three routes by which it is assembled.
- **The coloc arm regenerates `v2`, not the `v3` used by the published table.**
  `v3` is `v2` plus coloc recomputed for the 34 traits whose GWAS carried
  zero-SE or duplicate-ID SNPs (11,930,730 → 15,434,261 rows). That
  recomputation was run out of band and has no rule here.
- **The relaxed-threshold SMR run** that produces
  `tenk10k_phase1_sensitivity.*.parquet.gz` — the input to the sensitivity MR
  arm and to the preprocessing script — has no rule in this repository.
- **MR-link-2 is not vendored** and the Singularity image has no build recipe,
  so the `all_mrlink2` target is not currently runnable. IVW-LD, by contrast, is
  implemented in full in `sensitivity/rules/snakescripts/ivw_mr/run_ivw_corr.R`.
- **The BESD files** consumed by the SMR arm are staged, not built; the
  SAIGE-QTL → BESD conversion lives in the external `TenK10K_SMR` repository.

## Citation

If you use this code or data in your research, please consider citing

**Preprint**:

Henry, A., Senabouth, A., Tyebally, R., et al. Single-cell genetics identifies cell type-specific causal mechanisms in complex traits and diseases. medRxiv 2025.08.28.25334614 (2025) doi:10.1101/2025.08.28.25334614

[medRxiv link](https://www.medrxiv.org/content/10.1101/2025.08.28.25334614v2)

See [`CITATION.md`](CITATION.md) for citations on specific methods / datasets used in this study.


## Acknowledgments

- TenK10K study members and contributors
- Data providers and consortia (eQTLGen, IBDverse, GWAS consortia, STRING, Open Targets)
- Software developers & maintainers (Snakemake, R/Bioconductor, R Tidyverse, Python scientific computing ecosystem, MAGMA, scDRS, PLINK)
- Computational resources and support teams (Australia National Computing Infrastructure, Garvan Institute Data Science Platform)
