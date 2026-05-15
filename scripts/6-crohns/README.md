# 5-crohns: Crohn's Disease Case Study

This directory contains scripts for the Crohn's disease case study, which provides a detailed example of applying the causal inference framework to a specific disease with matched single-cell data.

## Overview

This section presents a comprehensive case study of Crohn's disease, integrating genetic associations, single-cell differential expression analysis, and causal inference. The analysis demonstrates how the TenK10K causal inference framework can be applied to understand disease mechanisms at cellular resolution.

## Contents

### Subdirectories

#### `deg/` - Differential Expression Gene Analysis
- `1-prepare_crohns_deg_pre-processed.R` - Prepare Crohn's DEG data for analysis
- `2-find_deg_mr_overlap_genes_both_colon_and_ti.R` - Find overlapping genes between DEG and MR results (colon and terminal ileum)
- `3-prepare_data_for_MR_comparison.py` - Prepare data for MR comparison analysis
- `4-plot_UMAP_cell_annotation.py` - Create UMAP plots with cell type annotations
- `5-plot_UMAP_genes_of_interest_counts.py` - Plot expression of genes of interest on UMAP

#### `figures/` - Figure Generation
- `Crohns_example_locus_zoom.R` - Create locus zoom plots for Crohn's disease examples
- `Figure5-combined_Crohns_figure.R` - Generate combined Figure 5 for manuscript
- `annotated_heatmap_only_canonicalanddeg.R` - Create annotated heatmaps for canonical and DEG results
- `barplot_comparison_with_eqtlgen_and_magma.R` - Comparison bar plots with eQTLGen and MAGMA
- `barplot_gene_numbers_by_celltype.R` - Bar plots showing gene numbers by cell type
- `crohns_deg_mr_comparison.R` - Comprehensive comparison of DEG and MR results

#### `prepare_data/` - Data Preparation
Scripts for preparing Crohn's disease-specific datasets

#### `supplementary/` - Supplementary Analyses
Additional analyses and supplementary figures for Crohn's disease case study
