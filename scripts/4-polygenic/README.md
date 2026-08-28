# 4-polygenic: Polygenic Enrichment Analysis

This directory contains scripts for polygenic enrichment analyses using single-cell Disease Relevance Score (scDRS) and scDeepID methods to identify cell-type-specific disease associations.

## Contents

### Subdirectories

#### `scDRS/`
Contains scripts for single-cell Disease Relevance Score analysis:
- `scdrs_main_supp.R` - Main scDRS analysis and figures for manuscript
- `scdrs_supp_sampling.R` - Supplementary sampling analyses for scDRS
- `compare_scdrs_aifi.R` - Replication against the AIFI Immune Health Atlas
- `aifi-labelTransfer.py` - Cell-type label transfer onto the AIFI atlas
- `README.md` - Detailed documentation, including which scripts a third party cannot run

#### `scDeepID/`
Documentation only. The scDeepID integration analysis was run from a separate
project tree with its own private module and input paths, so the script has been
removed rather than shipped in a form nobody could run. See `README.md`.

