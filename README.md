# TenK10K Phase1 effector genes identification manuscript

![Static Badge](https://img.shields.io/badge/version-preprint_v4-red)
[![medRxiv](https://img.shields.io/badge/medRxiv-10.1101%2F2025.08.28.25334614-B31B1B)](https://www.medrxiv.org/content/10.1101/2025.08.28.25334614v4)



## Study design
![](figures/biorender/study_design.png)

## Overview

This repository contains the code behind the manuscript *Single-cell genetics
identifies cell-type-specific effector genes across complex traits and
diseases*. It holds code and small metadata tables only, with no data or
results. Its two parts are:

- **Analysis pipelines** (`workflow/`, `sensitivity/`). These are Snakemake pipelines that turn GWAS and single-cell eQTL summary statistics into Mendelian randomisation (MR), sensitivity, colocalisation, MAGMA, scDRS, enrichment and Open Targets results. They need the TenK10K data and an HPC cluster.
- **Display items** (`scripts/`). These are the R, Python and Typst scripts that turn those results into every figure and table in the manuscript. Each item can be rebuilt with one command from the public summary statistics.

## Contents

| Path | Contents |
|---|---|
| `scripts/0-preprocess/` | Shared loader for the summary statistics, and helpers to fetch them |
| `scripts/1-overview/` … `scripts/8-crohns/` | One folder per Results subsection (see the table below); each has a README listing its scripts |
| `scripts/util/` | `write_table()`, which writes the supplementary tables |
| `scripts/Snakefile`, `Makefile` | Link every display item to the scripts and inputs that produce it |
| `workflow/` | Main analysis pipeline, plus the conda environments in `workflow/envs/` |
| `sensitivity/` | MR sensitivity pipeline ([IVW-MR](https://mrcieu.github.io/TwoSampleMR/index.html), [MR-link-2](https://github.com/adriaan-vd-graaf/mrlink2)) |
| `metadata/` | Cell types, traits, trait categories, gene annotation, supplementary-table column labels |
| `config/` | Settings for the display items (`figures.yaml`, `figure_inputs.txt`) and pipeline file maps |
| `figures/biorender/` | BioRender artwork for Fig. 1 |

## Reproducing the display items

**Requirements:** [conda](https://github.com/conda-forge/miniforge) and [Snakemake](https://snakemake.readthedocs.io) ≥ 8. Everything else is installed automatically.

```bash
git clone https://github.com/powellgenomicslab/tenk10k-effector-genes.git
cd tenk10k-effector-genes

make data        # download the summary statistics (once)
make fig-m2      # rebuild Figure 2  -> figures/fig-m2.pdf, figures/fig-m2.png
make all         # rebuild every figure and table
```

The first build installs a pinned R 4.4.1 environment, which takes a few minutes.

> [!NOTE]
> The summary-statistics download will be released with the published
> manuscript. Until then, `make data` reports that it is not yet available.

**Targets:**

| Command | Builds |
|---|---|
| `make fig-mN` | Main Figure N → `figures/fig-mN.pdf`, `.png` |
| `make fig-eN` | Extended Data Figure N → `figures/fig-eN.pdf`, `.png` |
| `make fig-sN` | Supplementary Figure N → `figures/fig-sN.pdf`, `.png` |
| `make tab-sN` | Supplementary Table N → `tables/tab-sN.tsv` |
| `make figs-m`, `figs-e`, `figs-s`, `tabs-s` | All items of one kind |
| `make list` | Every target, with its caption |

Dependencies between scripts are resolved automatically. For example, `make fig-m2` first builds the panels that Fig. 2 is assembled from. To build several items in parallel, add `CORES=4`.

### Display items by Results subsection

| Results subsection | Display items | Folder |
|---|---|---|
| Landscape of causal associations between gene expression and complex traits across immune cells | Fig. 1; Supp. Figs. 1–6; Supp. Tables 1–3 | [`1-overview`](scripts/1-overview/README.md) |
| Trait-level insights from cell-type-specific effector gene mapping | Fig. 2a–b; ED Fig. 1; Supp. Fig. 7; Supp. Tables 4–5 | [`2-trait`](scripts/2-trait/README.md) |
| Assessing horizontal pleiotropy via sensitivity analyses and colocalisation | Fig. 2c; ED Figs. 2–3; Supp. Fig. 8 | [`3-pleiotropy`](scripts/3-pleiotropy/README.md) |
| Single-cell eQTL enhances identification of effector genes | Fig. 3; ED Fig. 4; Supp. Figs. 9–11; Supp. Tables 6–7 | [`4-comparison`](scripts/4-comparison/README.md) |
| Identifying disease-relevant biological pathways | ED Fig. 5; Supp. Figs. 12–13; Supp. Table 8 | [`5-pathways`](scripts/5-pathways/README.md) |
| Polygenic enrichment reveals trait-relevant immune cell types | Fig. 4; ED Fig. 6; Supp. Figs. 14, 29; Supp. Tables 9–10 | [`6-polygenic`](scripts/6-polygenic/README.md) |
| Single-cell MR supports drug target discovery | Fig. 5; ED Fig. 8; Supp. Fig. 16; Supp. Tables 11–15 | [`7-drug`](scripts/7-drug/README.md) |
| Triangulation of evidence prioritises effector genes for Crohn's disease | Fig. 6a–c; ED Fig. 9; Supp. Figs. 17, 25, 30; Supp. Tables 16–21 | [`8-crohns`](scripts/8-crohns/README.md) |

Some items are not produced by code in this repository:
- **Fig. 1a** is BioRender artwork, provided in `figures/biorender/`.
- **Fig. 6** is assembled by hand from panels a–c (`make fig-m6`) and panel d.
- **Fig. 6d and Supp. Fig. 24** come from the [`TenK10K_causal_inference_vs_IBDverse`](TenK10K_causal_inference_vs_IBDverse) submodule.
- **ED Fig. 7 and Supp. Fig. 15** come from scDeepID, in a separate repository.

### Options

- **Container.** `make fig-m2 SDM="conda apptainer"` runs inside the published container image instead of building the environment locally. Apptainer is available on most HPC systems.
- **Your own software.** `make fig-m2 SDM=` uses the R and Python already on your `PATH`. The packages needed are listed in `workflow/envs/figures.yaml`.
- **Existing results tree.** If you already have the results (for example on NCI Gadi), run `make link-data DATA_DIR=/path/to/results/root` instead of `make data`. It links each input file listed in `config/figure_inputs.txt` and does not modify the source.
- **Raw-data steps.** Two kinds of step read individual-level data, which is not shared:
  - the LD calculations behind the GPX1 and JAZF1 locus plots;
  - the Crohn's colon UMAPs (Supp. Fig. 30).

  Their summary-level outputs are included in the download, so every other item builds without them. Where the raw data are available, add `ARGS="--config restricted=true"`.

### Regenerating the summary statistics

The pipelines in `workflow/` and `sensitivity/` need the TenK10K individual-level data and run on NCI Gadi under PBS; see [`workflow/README.md`](workflow/README.md) and [`sensitivity/README.md`](sensitivity/README.md). They are provided so that every step from raw data to summary statistics is documented. To check that a pipeline parses off the cluster:

```bash
snakemake --snakefile workflow/snakefile --workflow-profile none --list-rules
```

## Data Availability

The summary statistics needed to rebuild the display items will be deposited on
Zenodo / Hugging Face following manuscript publication; `make data` downloads
them. The sc-eQTL summary statistics are released separately as described in the main [TenK10K phase 1 sc-eQTL mapping study (Cuomo et al.)](https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2).

Not shared:
- Individual-level single-cell expression and genotypes of the TenK10K donors

Third-party data (please refer to the original publications for access):
- [Kong et al. 2023 single-cell transcriptomics from colon tissues](https://pubmed.ncbi.nlm.nih.gov/36720220/)
- [AIFI Immune Health Atlas (Gong et al. 2025)](https://doi.org/10.1038/s41586-025-09686-5)
- [IBDverse single-cell transcriptomics atlas (Alegbe T, Harris BT, et al.)](https://doi.org/10.1038/s41586-025-09686-5)

## Citation

If you use this code or data in your research, please consider citing

**Preprint**:

Henry, A., Senabouth, A., Tyebally, R., et al. Single-cell genetics identifies cell type-specific causal mechanisms in complex traits and diseases. medRxiv 2025.08.28.25334614 (2025) doi:10.1101/2025.08.28.25334614

[medRxiv link](https://www.medrxiv.org/content/10.1101/2025.08.28.25334614v2)

See *References* section in the manuscript for additional citations of methods and data sources.

## Acknowledgments

- TenK10K study members and contributors
- Data providers and consortia (eQTLGen, IBDverse, GWAS consortia, STRING, Open Targets)
- Software developers & maintainers (Snakemake, R/Bioconductor, R Tidyverse, Python scientific computing ecosystem, MAGMA, scDRS, PLINK)
- Computational resources and support teams (Australia National Computing Infrastructure, Garvan Institute Data Science Platform)
