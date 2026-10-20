# Snakemake Workflow for Tenk10k Causal Paper

## Overview
This Snakemake workflow processes and analyzes data for the Tenk10k causal inference paper, ensuring reproducibility and scalability.

## Features
- Data preprocessing
- Liftover of GWAS summary statistics
- Gene set enrichment analysis (gget Enrichr)
- SMR analysis
- Integration with Open Targets Platform
- MAGMA-based gene-level association analysis
- scDRS computation

## Requirements

- **Snakemake ≥ 8.** Only `profiles/default/config.v8+.yaml` is provided, and it
  uses the `cluster-generic` executor plugin, which is Snakemake 8+ syntax.
- Conda, for the environments in [`envs/`](envs/). All rules declare one, so
  pass `--use-conda`; the default profile already sets it.
- R ≥ 4.4 and Python ≥ 3.11 come from those environments.
- External tools **not** provided by the conda environments, which must be on
  `PATH`: `smr`, `magma`, `ldak`, `plink` (1.9), `tabix`/`htslib`, and `CrossMap`
  for the liftover chain. The multivariant coloc rules additionally expect NCI
  Gadi's `nci-parallel` module and PBS environment variables.

## Prerequisites you must supply

The workflow will parse without these, but rules will not run:

1. **`config/path/nci.yaml`** (or a host-specific equivalent) — a file map read at
   parse time by `rule init_file_nci`. A template with placeholder paths ships in
   the repository; replace the `source` entries with your own.
2. **`resources/ma/{trait}.ma`** — the curated GWAS inputs. This is an input, not
   an output; see the provenance note at the top of `rules/format_gwas.smk` and
   [`../scripts/0-preprocess/gwas/README.md`](../scripts/0-preprocess/gwas/README.md).
3. **Pre-built BESD files** for the sc-eQTL data, plus the GENCODE v44 GTF. Point
   `TENK10K_BESD_DIR` and `GENCODE_GTF` at them; see
   [`../scripts/0-preprocess/tenk10k-eqtl/README.md`](../scripts/0-preprocess/tenk10k-eqtl/README.md).
4. **Genotypes** under `resources/genotypes/` for LD, SMR and coloc.

## Usage

Run from the **working-directory root**, not from inside `workflow/`: every rule
declares its script paths relative to the root (for example
`workflow/rules/snakescripts/...`).

```bash
# Dry run - the quickest check that your configuration is complete
snakemake --snakefile workflow/snakefile \
          --profile workflow/profiles/default --dry-run

# A specific target
snakemake --snakefile workflow/snakefile \
          --profile workflow/profiles/default \
          results/aggregate/tenk10k_phase1.msmr.parquet.gz

# Locally, with N cores instead of submitting to PBS
snakemake --snakefile workflow/snakefile --use-conda --cores 8 <target>
```

The default profile submits to PBS on NCI Gadi (project `fy54`, `cluster-generic`).
Edit `profiles/default/config.v8+.yaml` — in particular `project`, `storage` and
`email` — or supply your own profile.

## Workflow Steps

### 1. GWAS Formatting
- **Scripts:** `rules/snakescripts/format_gwas/*.R`
- **Purpose:** Standardize GWAS summary statistics, liftover coordinates, harmonize alleles, filter variants.
- **Output:** Formatted files in `resources/pipeline_ma/` (a staging area; see the provenance note in `rules/format_gwas.smk` for how these reach `resources/ma/`)

### 2. eQTL Preparation
- **Scripts:** `rules/snakescripts/prep_besd_chr/*.sh`, `prep_smr_input/*.R`
- **Purpose:** Prepare BESD files, probe lists, and p-value thresholds for SMR.
- **Output:** `resources/besd/`, `resources/smr/`

### 3. SMR Analysis
- **Scripts:** `rules/snakescripts/run_smr.sh`, `smr_locus.sh`
- **Purpose:** Run SMR per chromosome/cell type/phenotype, extract locus information.
- **Output:** `results/smr/`, `results/smr_locus/`

### 4. MAGMA Gene-Level Association
- **Scripts:** `rules/snakescripts/aggregate/magma.R`, `magma_format_output.R`
- **Purpose:** Aggregate MAGMA results, annotate genes, FDR correction.
- **Output:** `results/magma/`

### 5. Gene Set Enrichment (Enrichr)
- **Scripts:** `rules/snakescripts/enrichment/gget_enrichr.py`, `gget_enrichr_pheno.py`
- **Purpose:** Perform gene set enrichment using Enrichr via gget.
- **Output:** `results/enrichment/`

### 6. scDRS Computation
- **Scripts:** `rules/snakescripts/scdrs/*`
- **Purpose:** Prepare covariates, regress out confounders, compute scDRS scores.
- **Output:** `resources/scdrs/`, `results/scdrs/`
