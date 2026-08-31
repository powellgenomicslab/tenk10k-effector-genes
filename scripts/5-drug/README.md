# 5-drug: Drug Target Support Analysis

This directory contains scripts for analyzing existing evidence of clinical development for causal genes identified through Mendelian Randomization analyses using data from [Open Targets Platform](https://platform.opentargets.org/).

## Overview

This section performs drug target enrichment analyses to identify whether genes with causal effects on disease traits are enriched for known drug targets. 

## Contents

### Scripts
- `otp_combined.R` - **canonical** drug enrichment and target-support analysis,
  against Open Targets Platform release **26.03**. It:
  - analyses enrichment of drug targets among causal genes
  - compares categories of target support and clinical phase
  - produces the drug-target figures and supplementary tables
- `otp_power_projection.R` - projects the statistical power of the
  target-concordance analysis, also against release 26.03

An earlier version of this analysis ran against Open Targets release 25.06 and
emitted the same four supplementary tables. It has been removed to leave one
canonical version; see the git history if you need it.

### Prerequisites not in this repository
- `metadata/otp.xlsx` - therapeutic-area and clinical-phase mappings
- `resources/misc/chembl_action_type.tsv`
- `results/otp/26.03/otp_concordance_by_phase.tsv` - written by a commented-out
  block in `otp_combined.R`

