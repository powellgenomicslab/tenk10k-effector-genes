# TenK10K Phase1 effector genes identification manuscript

![Static Badge](https://img.shields.io/badge/version-preprint_v4-red)
[![medRxiv](https://img.shields.io/badge/medRxiv-10.1101%2F2025.08.28.25334614-B31B1B)](https://www.medrxiv.org/content/10.1101/2025.08.28.25334614v4)



## Study design
![](figures/biorender/study_design.png)

## Overview

This repository holds the code for the TenK10K effector genes manuscript: the
Snakemake pipelines that turn raw data into summary statistics, and the scripts
that turn those summary statistics into every figure and table in the paper.
It contains code and small metadata only; no data or results.

## Directory structure

- **`scripts/`**: figure and table scripts, one folder per Results subsection
  of the manuscript.

  | Folder | Results subsection |
  |---|---|
  | [`0-preprocess/`](scripts/0-preprocess/README.md) | Shared loader for the summary statistics, plus helpers for the input bundle |
  | [`1-overview/`](scripts/1-overview/README.md) | Landscape of causal associations between gene expression and complex traits across immune cells |
  | [`2-trait/`](scripts/2-trait/README.md) | Trait-level insights from cell-type-specific effector gene mapping |
  | [`3-pleiotropy/`](scripts/3-pleiotropy/README.md) | Assessing horizontal pleiotropy in single-cell eQTL Mendelian randomisation via sensitivity analyses and colocalisation |
  | [`4-comparison/`](scripts/4-comparison/README.md) | Single-cell eQTL enhances identification of effector genes for complex traits |
  | [`5-pathways/`](scripts/5-pathways/README.md) | Identifying disease-relevant biological pathways using single-cell genetics |
  | [`6-polygenic/`](scripts/6-polygenic/README.md) | Polygenic enrichment reveals trait-relevant immune cell types and gene programs |
  | [`7-drug/`](scripts/7-drug/README.md) | Single-cell MR supports drug target discovery and clinical translation |
  | [`8-crohns/`](scripts/8-crohns/README.md) | Single-cell genetics coupled with triangulation of evidence prioritises effector genes for Crohn's disease |
  | `util/` | `write_table()`, which writes supplementary tables |

  `scripts/Snakefile` connects every script to the manuscript item it produces.
- **`workflow/`**: the main Snakemake pipeline, from GWAS and sc-eQTL summary statistics to SMR, MAGMA, colocalisation, scDRS, enrichment and Open Targets results. `workflow/envs/` holds the conda environments, including `figures.yaml` for `scripts/`.
- **`sensitivity/`**: the Snakemake pipeline for the MR sensitivity analyses, [IVW-MR](https://mrcieu.github.io/TwoSampleMR/index.html) and [MR-link-2](https://github.com/adriaan-vd-graaf/mrlink2).
- **`metadata/`**, **`config/`**: small tracked tables and settings, including cell types, traits, gene annotation, supplementary-table column labels, file maps and `figures.yaml`.
- **`figures/biorender/`**: BioRender artwork used in Fig. 1. Generated figures are written to `figures/` and are not tracked.

## Reproducing the figures and tables

Each manuscript item has one command:

| Command | Builds |
|---|---|
| `make fig-m3` | Main Figure 3 → `figures/fig-m3.pdf`, `.png` |
| `make fig-e3` | Extended Data Figure 3 → `figures/fig-e3.pdf`, `.png` |
| `make fig-s3` | Supplementary Figure 3 → `figures/fig-s3.pdf`, `.png` |
| `make tab-s3` | Supplementary Table 3 → `tables/tab-s3.tsv` |
| `make figs-m`, `figs-e`, `figs-s`, `tabs-s`, `all` | groups of items |
| `make list` | every target with its caption |

The Makefile is a thin wrapper around [`scripts/Snakefile`](scripts/Snakefile). Snakemake resolves the dependencies between scripts. For example, Fig. 2 is assembled from panels made by three other scripts.

**1. Software.** You need Snakemake ≥ 8 and conda (Miniforge is recommended).
- `make fig-m3` builds the pinned `figures` environment on first use, with exact package versions from `workflow/envs/figures.<platform>.pin.txt`.
- Once the container image is published, `make fig-m3 SDM="conda apptainer"` runs inside it instead. Use this on Gadi.
- `SDM=` uses whatever R and Python are already on your `PATH`.

**2. Inputs.** The scripts read summary statistics only. Use one of:
- `make data` downloads the public bundle (Zenodo / Hugging Face; the link is in `config/figures.yaml`) and checks every file against its checksum.
- `make link-data DATA_DIR=/g/data/fy54/analysis/tenk10k-causal` links the same files from an existing results tree. It links file by file, so anything the scripts write stays in this checkout.

`config/figure_inputs.txt` lists every input, and `make manifest` regenerates it from the Snakefile.

**3. Build.** `make all` builds every item that can be made from summary statistics.

A few steps read raw single-cell data or individual-level genotypes:
- the Crohn's colon UMAPs (Supplementary Fig. 30);
- the LD calculations behind the GPX1 and JAZF1 locus objects.

These steps run only with `make <target> ARGS="--config restricted=true"`, where the raw data are available. Otherwise their summary-level outputs come from the bundle.

Off Gadi, the main pipeline needs `--workflow-profile none`, because its default profile submits jobs to the PBS cluster:
```bash
snakemake --snakefile workflow/snakefile --workflow-profile none --list-rules
```

### Manuscript items and their scripts

Folder numbers match the Results subsections.

| Item | Target | Folder |
|---|---|---|
| Fig. 1 | `fig-m1` | 1 (BioRender panel a and Typst layout) |
| Fig. 2 | `fig-m2` | 2 (panels a–b) and 3 (panel c) |
| Fig. 3 | `fig-m3` | 4 |
| Fig. 4 | `fig-m4` | 6 |
| Fig. 5 | `fig-m5` | 7 |
| Fig. 6a–c | `fig-m6` | 8 (panels only; the figure is assembled by hand) |
| ED Figs. 1, 2–3, 4, 5, 6, 8, 9 | `fig-e1` … `fig-e9` | 2, 3, 4, 5, 6, 7, 8 |
| Supp. Figs. 1–6 | `fig-s1` … `fig-s6` | 1 |
| Supp. Fig. 7 | `fig-s7` | 2 |
| Supp. Fig. 8 | `fig-s8` | 3 |
| Supp. Figs. 9–11 | `fig-s9` … `fig-s11` | 1 (Fig. 9) and 4 |
| Supp. Figs. 12–13 | `fig-s12`, `fig-s13` | 5 |
| Supp. Figs. 14, 29 | `fig-s14`, `fig-s29` | 6 |
| Supp. Fig. 16 | `fig-s16` | 7 |
| Supp. Figs. 17, 25, 30 | `fig-s17`, `fig-s25`, `fig-s30` | 8 |
| Supp. Tables 1–21 | `tab-s1` … `tab-s21` | 1–8 |

Items with no script in this repository:
- **Fig. 1a**: BioRender artwork, shipped in `figures/biorender/`.
- **Fig. 6d and Supp. Fig. 24**: made in the [`TenK10K_causal_inference_vs_IBDverse`](TenK10K_causal_inference_vs_IBDverse) submodule.
- **ED Fig. 7 and Supp. Fig. 15**: made with scDeepID, in a separate repository.
- **Supp. Figs. 18–23 and 26–28**: not yet mapped to a script.

## Data Availability

Aggregated summary statistics — the MR, colocalisation and related result
tables — will be deposited on Hugging Face following manuscript
publication. The sc-eQTL summary statistics are released separately as described in the main [TenK10K phase 1 sc-eQTL mapping study (Cuomo et al.)](https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2).

Data that are not currently included in this repository:
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
