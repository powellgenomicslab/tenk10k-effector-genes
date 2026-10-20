# TenK10K Phase1 effector genes identification manuscript

[![medRxiv](https://img.shields.io/badge/medRxiv-10.1101%2F2025.08.28.25334614-B31B1B)](https://doi.org/10.1101/2025.08.28.25334614)
[![Data](https://img.shields.io/badge/data-10.5281%2Fzenodo.23406585-1682D4)](https://doi.org/10.5281/zenodo.23406585)



## Study design
![](figures/biorender/study_design.png)

## Overview

This repository contains the code behind the manuscript *Single-cell genetics
identifies cell-type-specific effector genes across complex traits and
diseases*. It holds code and small metadata tables only, with no data or
results. Its two parts are:

- **Analysis pipelines** (`workflow/`, `sensitivity/`). These are Snakemake pipelines that turn GWAS and single-cell eQTL summary statistics into Mendelian randomisation (MR), sensitivity, colocalisation, MAGMA, scDRS, enrichment and Open Targets results. They need the TenK10K data (not included due to size).
- **Display items** (`scripts/`). These are the R, Python and Typst scripts that turn those results into every figure and table in the manuscript. Each item has a single command that rebuilds it from summary-level inputs (released separately on Zenodo, doi:[10.5281/zenodo.23406585](https://doi.org/10.5281/zenodo.23406585); see [Data Availability](#data-availability)).

## Main results

> [!IMPORTANT]
> The study's main results are one file in the Zenodo data record
> ([doi:10.5281/zenodo.23406585](https://doi.org/10.5281/zenodo.23406585)):
> **`tenk10k_phase1_effector_genes.v5.parquet`**.

It holds the single-cell eQTL Mendelian randomisation results for every gene,
cell type and trait tested in TenK10K Phase 1 (7.5 million rows): multi-SNP SMR
with local false discovery rates, the HEIDI, IVW-LD and MR-link-2 sensitivity
analyses, coloc and SuSiE-based colocalisation, and the evidence tiers used in
the manuscript (`mr`, `mr_sens`, `mr_coloc`, `mr_sens_coloc`). Each column is
described in `tenk10k_phase1_effector_genes.v5.md`, next to it in the record.
The same descriptions are stored in the file's own metadata, so it documents
itself.

```r
df <- arrow::read_parquet("tenk10k_phase1_effector_genes.v5.parquet")
subset(df, mr_sens_coloc)   # MR, sensitivity and colocalisation all supported
```

```python
import pandas as pd
df = pd.read_parquet("tenk10k_phase1_effector_genes.v5.parquet")
df[df.mr_sens_coloc]
```

The file is made from `results/preprocessed/tenk10k_phase1.v5.parquet.gz` by
the rule `main_results_release` (`workflow/rules/misc.smk`): it has every
column of that file, renamed and documented in
[`workflow/config/main_results_columns.tsv`](workflow/config/main_results_columns.tsv).

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

### Quick start

You need [conda](https://github.com/conda-forge/miniforge), and about 32 GB
of memory and 10 GB of disk space.

```bash
# 1. Install Snakemake (once); activate it whenever you run make
conda create -n snakemake -c conda-forge -c bioconda 'snakemake>=9'
conda activate snakemake

# 2. Get the code and the data (once)
git clone https://github.com/powellgenomicslab/tenk10k-effector-genes.git
cd tenk10k-effector-genes
make data

# 3. Build a figure -> figures/fig-m2.pdf, figures/fig-m2.png
make fig-m2
```

`make all` builds every figure and table, and `make list` shows each target
with its caption. To run several scripts at once, add `CORES=N`. Each script
can use up to about 16 GB of memory, so choose N no larger than your memory in
GB divided by 16 (for example `CORES=4` with 64 GB). Rebuilding everything takes about
an hour of computing: from a fresh clone on NCI Gadi, `make data` took under a
minute and `make all` with `CORES=4` took 21 minutes and at most 47 GB. The first
build installs a pinned conda environment, which takes several minutes; later
builds reuse it.

> [!NOTE]
> The data record ([doi:10.5281/zenodo.23406585](https://doi.org/10.5281/zenodo.23406585)) is published, but its
> files stay restricted until the manuscript is published. Until then, `make data`
> cannot download them.

### Using the container instead of conda

If [Apptainer](https://apptainer.org) (formerly Singularity) is available, as
on most Linux workstations and HPC systems, you can run the scripts inside a
ready-made container image instead of installing the conda environment. Add
`SDM="conda apptainer"` to any build:

```bash
make fig-m2 SDM="conda apptainer"
```

The first such build downloads the image (1.5 GB); later builds reuse it.

### Running on an HPC cluster

On most clusters, compute nodes have no internet access and login nodes are
not meant for heavy work. Split the work in two:

1. **Download**, once, wherever your site allows jobs with internet access
   (for example a data-transfer queue). Load Apptainer first if it is a
   module (`module load apptainer`), and allow 16 GB of memory and about 6 GB of
   temporary disk (`make data` unpacks the archive in `$TMPDIR`, which on some
   clusters is a small per-job disk; on NCI Gadi, request `-l jobfs=20GB`):

   ```bash
   make data
   make pull          # container image only; runs no scripts
   ```

2. **Build** in a batch or interactive job on a compute node, with about
   16 GB of memory per core (here 4 cores and 64 GB):

   ```bash
   conda activate snakemake
   make all SDM="conda apptainer" CORES=4
   ```

If you set `APPTAINER_CACHEDIR`, use the same value in both steps, so that the
build finds the downloaded image. The container can see the repository
directory; if your data are linked from elsewhere (see `make link-data` under
[Options](#options)), make that path visible with
`ARGS='--apptainer-args "-B /path/to/data"'`.

**Targets:**

| Command | Builds |
|---|---|
| `make fig-mN` | Main Figure N → `figures/fig-mN.pdf`, `.png` |
| `make fig-eN` | Extended Data Figure N → `figures/fig-eN.pdf`, `.png` |
| `make fig-sN` | Supplementary Figure N → `figures/fig-sN.pdf`, `.png` |
| `make tab-sN` | Supplementary Table N → `tables/tab-sN.tsv` |
| `make figs-m`, `figs-e`, `figs-s`, `tabs-s` | All items of one kind |
| `make list` | Every target, with its caption |
| `make figures/panels/<folder>/<file>` | One output file by path, such as a panel saved on its own (e.g. `make figures/panels/1-overview/panel_b.svg`). Scripts that draw several panels write them all |

Dependencies between scripts are resolved automatically. For example, `make fig-m2` first builds the panels that Fig. 2 is assembled from.

### Display items by Results subsection

| Results subsection | Display items | Folder |
|---|---|---|
| Landscape of causal associations between gene expression and complex traits across immune cells | Fig. 1; Supp. Figs. 1–6; Supp. Tables 1–3 | [`1-overview`](scripts/1-overview/README.md) |
| Trait-level insights from cell-type-specific effector gene mapping | Fig. 2a–b; ED Fig. 1; Supp. Fig. 7; Supp. Tables 4–5 | [`2-trait`](scripts/2-trait/README.md) |
| Assessing horizontal pleiotropy via sensitivity analyses and colocalisation | Fig. 2c; ED Figs. 2–3; Supp. Fig. 8 | [`3-pleiotropy`](scripts/3-pleiotropy/README.md) |
| Single-cell eQTL enhances identification of effector genes | Fig. 3; ED Fig. 4; Supp. Figs. 9–11; Supp. Tables 6–7 | [`4-comparison`](scripts/4-comparison/README.md) |
| Identifying disease-relevant biological pathways | ED Fig. 5; Supp. Figs. 12–13; Supp. Table 8 | [`5-pathways`](scripts/5-pathways/README.md) |
| Polygenic enrichment reveals trait-relevant immune cell types | Fig. 4; ED Figs. 6–7; Supp. Figs. 14–15, 21–24; Supp. Tables 9–10 | [`6-polygenic`](scripts/6-polygenic/README.md) |
| Single-cell MR supports drug target discovery | Fig. 5; ED Fig. 8; Supp. Fig. 16; Supp. Tables 11–15 | [`7-drug`](scripts/7-drug/README.md) |
| Triangulation of evidence prioritises effector genes for Crohn's disease | Fig. 6; ED Figs. 9–10; Supp. Figs. 17, 25 (raw data only); Supp. Tables 16–21 | [`8-crohns`](scripts/8-crohns/README.md) |

Some items are not produced by code in this repository:
- **Fig. 1a** is BioRender artwork, provided in `figures/biorender/`.
- **Fig. 6** is assembled by `make fig-m6` (`scripts/8-crohns/fig6.typ`); panel d is ported from the IBDverse comparison (`scripts/8-crohns/plot_figures/ibdverse_blood_enrichment.R`).
- **Supp. Figs. 18–19** come from the [`TenK10K_causal_inference_vs_IBDverse`](TenK10K_causal_inference_vs_IBDverse) submodule.
- **Supp. Fig. 20** (genes tested per cell type) is not yet linked to a script.

### Options

- **Your own software.** `make fig-m2 SDM=` uses the R and Python already on your `PATH`. The packages needed are listed in `workflow/envs/figures.yaml`.
- **Rebuilding the container image.** The image is built from the `Dockerfile`, which `make dockerfile` regenerates from the conda environments. Publish each rebuild under a new tag and update `container:` in `config/figures.yaml`: Snakemake caches images under a name derived from the tag, so a rebuilt image pushed under an existing tag is ignored wherever the old one has been downloaded.
- **Existing results tree.** If you already have the results tree, run `make link-data DATA_DIR=/path/to/results/root` instead of `make data`. It links each input file listed in `config/figure_inputs.txt` and does not modify the source.
- **Not rebuildable from the download.** Supp. Fig. 25 (UMAPs of the Kong et al. 2023 colon immune cells) needs their single-cell data, which is not shared, so it is drawn by the analysis pipeline (`workflow/rules/crohns.smk`); `make fig-s25` says so. Everything else is built from the download alone: steps that start from raw or unshared data (genotype LD, the Crohn's literature and DEG tables, the AIFI label transfer, the scDeepID cell function scores) are in `workflow/`, and only their results are in the download.

### Regenerating the summary statistics

The pipelines in `workflow/` and `sensitivity/` need the TenK10K individual-level data and run on NCI Gadi under PBS; see [`workflow/README.md`](workflow/README.md) and [`sensitivity/README.md`](sensitivity/README.md). They are provided so that every step from raw data to summary statistics is documented. To check that a pipeline parses off the cluster:

```bash
snakemake --snakefile workflow/snakefile --workflow-profile none --list-rules
```

## Data Availability

The data needed to rebuild the display items are deposited on Zenodo as a
dataset with its own DOI, [10.5281/zenodo.23406585](https://doi.org/10.5281/zenodo.23406585), separate from this code's
release; its files are restricted until the manuscript is published. `make data` downloads the archive and checks every file
against its checksum. The archive's `INVENTORY.tsv` lists each file with its
source and licence. [`RELEASING.md`](RELEASING.md) describes how the code and
data releases are linked. The sc-eQTL summary statistics are released separately as described in the main [TenK10K phase 1 sc-eQTL mapping study (Cuomo et al.)](https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2).

Not shared:
- Individual-level single-cell expression and genotypes of the TenK10K donors

Third-party data (please refer to the original publications for access):
- [Kong et al. 2023 single-cell transcriptomics from colon tissues](https://pubmed.ncbi.nlm.nih.gov/36720220/)
- [AIFI Immune Health Atlas (Gong et al. 2025)](https://doi.org/10.1038/s41586-025-09686-5)
- [IBDverse single-cell transcriptomics atlas (Alegbe T, Harris BT, et al.)](https://doi.org/10.6019/S-BSST2922)

## Contributors

The code in this repository was written by:

- Albert Henry ([@alhenry](https://github.com/alhenry))
- Anne Senabouth ([@asenabouth](https://github.com/asenabouth))
- Blake Bowen ([@blake-bowen](https://github.com/blake-bowen))
- Jayden Fan ([@camwew](https://github.com/camwew))
- Oscar Dong ([@OscarADong](https://github.com/OscarADong))
- Peter Allen ([@PeterCAllen](https://github.com/PeterCAllen))
- Lawrence Huang ([@Hao-Lawrence-Huang](https://github.com/Hao-Lawrence-Huang))
- Rika Tyebally ([@riktye](https://github.com/riktye))
- Bradley T. Harris ([@BradleyH017](https://github.com/BradleyH017))

> [!NOTE]
> **Use of AI.** Claude (Anthropic) was used to organise and curate this
> repository for release: arranging the code by manuscript section, writing the
> build pipeline and documentation, and checking that every display item can be
> rebuilt from the released data. The analyses, and the code that implements
> them, are the work of the people named above.

## Citation

If you use this code or data in your research, please consider citing

**Preprint**:

Henry, A., Senabouth, A., Tyebally, R., et al. Single-cell genetics identifies cell type-specific causal mechanisms in complex traits and diseases. medRxiv 2025.08.28.25334614 (2025) doi:10.1101/2025.08.28.25334614

[medRxiv link](https://doi.org/10.1101/2025.08.28.25334614)

See *References* section in the manuscript for additional citations of methods and data sources.

## Acknowledgments

- TenK10K study members and contributors
- Data providers and consortia (eQTLGen, IBDverse, GWAS consortia, STRING, Open Targets)
- Software developers & maintainers (Snakemake, R/Bioconductor, R Tidyverse, Python scientific computing ecosystem, MAGMA, scDRS, PLINK)
- Computational resources and support teams (Australia National Computing Infrastructure, Garvan Institute Data Science Platform)
