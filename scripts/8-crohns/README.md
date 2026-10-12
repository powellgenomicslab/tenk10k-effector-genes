# 6-crohns: Crohn's Disease Case Study

This directory contains scripts for the Crohn's disease case study, which provides a detailed example of applying the causal inference framework to a specific disease with matched single-cell data.

## Overview

This section presents a comprehensive case study of Crohn's disease, integrating genetic associations, single-cell differential expression analysis, and causal inference. The analysis demonstrates how the TenK10K causal inference framework can be applied to understand disease mechanisms at cellular resolution.

## Contents

### `crohns_case_study/` — the canonical analysis

- `prepare_data/` - postprocessing of the MR results and the external-literature
  gene annotation table
- `deg/` - Kong et al. 2023 differential expression, its overlap with the MR
  results, and the UMAP panels
- `plot_figures/` - the figure panels (gene annotation, shared immune diseases,
  annotated heatmap, MR/DEG scatter, GPX1 and JAZF1 loci)
- `supplementary/` - the Crohn's supplementary workbook

### `figures/` and `supplementary/` — earlier generation

Retained for reference. `supplementary/1-prepare_supplementary_table_overall_results.R`
is still the sole producer of `crohns_summary_gene_count_df.RData`, which the
canonical supplementary script loads.

> [!NOTE]
> Earlier generations of `prepare_data/` and `deg/` have been removed so that
> there is one canonical version of each step; see the git history if needed.
> The locus-zoom scripts need donor-derived LD matrices and genotypes, which are
> not part of the public data release.
