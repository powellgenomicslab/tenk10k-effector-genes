# 2-mr: Mendelian Randomization Analysis

This directory contains scripts for Mendelian Randomization (MR) analyses and visualization of results

## Contents

### Scripts
- `mr_results_main.R` - Main figures and analyses for MR results:
  - Creates Figure 2 showing MR overview across cell types and traits
  - Generates tally plots comparing MR-only vs MR & GWAS associations
  - Compares TenK10K single-cell results with eQTLGen bulk results
  - Creates comprehensive visualization panels for manuscript

- `mr_results_supp.R` - Supplementary MR analyses:
  - Additional statistical comparisons
  - Extended visualizations for supplementary materials
  - Detailed breakdowns by cell type and trait categories

## Dependencies

### R packages
```r
library(patchwork)      # Combining multiple plots
library(ragg)           # High-quality graphics device
library(scales)         # Scale functions for visualization
library(paletteer)      # Color palettes
library(geomtextpath)   # Text along paths in plots
# Plus all dependencies from 0-preprocess/preprocess_results.R
```

## Input Files

Both scripts source `scripts/0-preprocess/preprocess_results.R`, which provides:
- Filtered mSMR results from TenK10K and eQTLGen
- Cell type and trait mappings
- Gene annotations and categories


## Color Schemes

The scripts use consistent color schemes:
- `mr_only` - Light blue (#A4ABB0FF) for MR-only associations
- `mr_other` - Dark blue (#4C6C94FF) for MR & GWAS associations
