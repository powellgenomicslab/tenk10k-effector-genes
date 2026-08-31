# 0-preprocess: Data Preprocessing

This directory contains scripts and data for preprocessing GWAS summary statistics, eQTL, and the main MAGMA & Mendelian randomization analysis results for downstream analyses.

## Contents

### Scripts
- `preprocess_results.R` - Main preprocessing script that:
  - Loads cell type and trait mappings from metadata
  - Processes MAGMA gene-level association results
  - Filters and harmonizes MR results from TenK10K and eQTLGen
  - Creates a unified gene universe for downstream analyses
  - Applies quality control filters (removes missing p-values, zero effect sizes)

The script applies several QC filters:
- Removes genes with missing SMR p-values
- Excludes associations with zero GWAS or SMR effect sizes
- Creates gene universe based on intersection of genes available in both MAGMA and TenK10K analyses
- Applies trait inclusion criteria from metadata

### Subdirectories

Both hold a README only. The code they describe lives in the Snakemake pipeline,
and each README points at it:
- [`gwas/`](gwas/README.md) - where GWAS harmonisation happens, and how
  `resources/ma/` is assembled
- [`tenk10k-eqtl/`](tenk10k-eqtl/README.md) - how SAIGE-QTL output enters the
  pipeline, and what you must supply

## Dependencies

### R packages
```r
library(data.table)    # Fast data manipulation
library(tidyverse)     # Data wrangling and visualization
library(arrow)         # Parquet file handling
library(fs)            # File system operations
library(qvalue)        # Multiple testing correction
```

## Input Files

Tracked in this repository (no action needed):
- `metadata/cell.tsv` - cell-type mapping, colours and ordering
- `metadata/trait.tsv` - trait definitions, sample sizes and inclusion flags
- `metadata/trait_category.tsv` - trait category ordering and colours
- `metadata/gencode.v44.gene_type.tsv.gz` - gene annotation

From the data release, under the working directory:
- `results/preprocessed/tenk10k_phase1.v5.parquet.gz` - the integrated MR /
  sensitivity / coloc evidence table (see [`../../CHANGELOG.v5.md`](../../CHANGELOG.v5.md))
- `results/aggregate/tenk10k_phase1.magma.gz.parquet` - MAGMA results
- `results/aggregate/eqtlgen2020.msmr.parquet.gz` - eQTLGen mSMR results
- `results/aggregate/coloc/tenk10k_phase1.coloc.v3.parquet.gz` and
  `tenk10k_phase1.mvcoloc.parquet.gz` - colocalisation results

## Usage

Run the preprocessing script from the repository root:

```bash
Rscript scripts/0-preprocess/preprocess_results.R
```

**Note:** Ensure that the Snakemake workflow has been executed first to generate the required input files in `results/aggregate/`.

