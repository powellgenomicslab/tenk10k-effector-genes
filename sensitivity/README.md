# Sensitivity Analysis Pipeline

This Snakemake pipeline runs two sensitivity MR methods — **IVW-LD** and **MR-link-2** — on significant eQTL–trait associations identified from the TenK10K Phase 1 SMR-multi analysis.

## Overview

The pipeline:
1. Reads SMR-multi results and generates per-trait phenotype lists (cell types with significant associations at lfdr < 0.05)
2. Runs IVW with LD correction (IVW-LD) for each trait × cell type combination
3. Runs MR-link-2 for each trait × cell type × gene combination
4. Aggregates results into compressed Parquet outputs

## Directory Structure

```
sensitivity/
├── Snakefile                          # Main workflow entry point
├── config/
│   └── gadi.yaml                      # NCI Gadi cluster configuration
├── rules/
│   ├── ivw_mr.smk                     # IVW-LD rules
│   ├── mrlink2.smk                    # MR-link-2 rules
│   └── snakescripts/
│       ├── ivw_mr/
│       │   ├── run_ivw_corr.R         # Run IVW-LD per trait × cell type
│       │   ├── process_ivw_corr.R     # Process per-trait IVW results
│       │   └── make_ivw_corr_parquet.R # Combine into final Parquet
│       └── mrlink2/
│           ├── preprocess_mrlink2_outcome.R        # Prepare GWAS outcome data
│           ├── preprocess_mrlink2_exposure.R       # Prepare eQTL exposure data
│           ├── preprocess_mrlink2_exposure_traits.R # Filter genes per trait
│           ├── process_mrlink2_results.R           # Aggregate per-trait results
│           └── combine_mrlink2_results.R           # Combine into final Parquet
└── scripts/
    └── preparePhenotypeList.R         # Checkpoint: generate per-trait cell type lists
```

## Requirements

- Snakemake
- Singularity (jobs run inside `mr.sif` container)
- NCI Gadi access (project `fy54`)
- Input data on `/g/data/fy54/` (see [Configuration](#configuration))
- Singularity image: [Docker Hub](https://hub.docker.com/repository/docker/asenabouth/mendelianrandomization)

## Configuration

The pipeline auto-detects whether it is running on a cluster (PBS or SLURM) or locally:

- **Cluster (Gadi):** loads `sensitivity/config/gadi.yaml`
- **Local:** loads `sensitivity/config/local.yaml`

Both are evaluated at parse time, so the relevant one must exist before any rule
can be considered. A documented `local.yaml` template now ships in the
repository; edit the paths in it rather than creating one from scratch.

> [!NOTE]
> `gadi.yaml` conflates two different things: Snakemake **config** keys
> (`workflow`, `key_files`, `container_paths`) and Snakemake **profile** settings
> (`cluster`, `default-resources`, `jobs`, `use-singularity`, `singularity-args`).
> The `Snakefile` loads it via `configfile:`. The `--profile` invocations below
> are reproduced as originally documented, but note that `--profile` expects a
> *directory*, and the `cluster:` key is the pre-Snakemake-8 interface — so these
> commands are unlikely to work as written on a current Snakemake. Prefer passing
> the file with `--configfile` and supplying execution settings separately.

Key paths in `gadi.yaml`:

| Key | Description |
|-----|-------------|
| `key_files.smr_results` | SMR-multi results Parquet (input) |
| `key_files.eqtl_dir` | eQTL summary statistics directory |
| `workflow.ivw_output_dir` | IVW-LD output directory (scratch) |
| `workflow.mrlink2_output_dir` | MR-link-2 output directory (scratch) |
| `mr_singularity_image` | Path to `mr.sif` Singularity container |

Singularity bind mounts (configured in `gadi.yaml`):

| Container path | Host path |
|---|---|
| `/data` | `/g/data/fy54/analysis/tenk10k-causal` |
| `/genotypes` | `/g/data/fy54/genotypes/common_variants/2024-12_Freeze/plink` |
| `/workspace` | `/scratch/fy54/as8574/analysis/tenk10k-causal` |
| `/eqtl` | `/g/data/fy54/results/eqtl/saige_qtl/december24_freeze` |

## Running the Pipeline

The pipeline must be run from the **repository root** (not from inside `sensitivity/`), because Snakemake rule scripts use paths relative to the root. `sensitivity/Snakefile` is the only entry point; pass it with `--snakefile`.

### Run IVW-LD only

```bash
snakemake --snakefile sensitivity/Snakefile \
    --profile sensitivity/config/gadi.yaml \
    all_ivw
```

### Run MR-link-2 only

```bash
snakemake --snakefile sensitivity/Snakefile \
    --profile sensitivity/config/gadi.yaml \
    all_mrlink2
```

### Dry run (check workflow without executing)

```bash
snakemake --snakefile sensitivity/Snakefile \
    --profile sensitivity/config/gadi.yaml \
    --dry-run all_ivw
```

### Submit to PBS queue (Gadi)

The `gadi.yaml` config includes PBS `cluster` settings. Pass it as a profile:

```bash
snakemake --snakefile sensitivity/Snakefile \
    --profile sensitivity/config/gadi.yaml \
    all_ivw all_mrlink2
```

## Workflow Steps

### Checkpoint: `prepare_phenotype_list`

Reads the SMR-multi results Parquet and filters to associations with lfdr < 0.05. Outputs:
- `sensitivity/metadata/phenotypes/{trait}.txt` — cell types with significant hits per trait
- `sensitivity/metadata/phenotype_list.txt` — list of all traits to process

This checkpoint gates all downstream rules; trait lists are resolved dynamically.

### IVW-LD (`all_ivw`)

| Rule | Description |
|------|-------------|
| `run_ivw` | Runs IVW-LD for each trait across its cell types (8 CPUs, 64 GB, express queue) |
| `process_ivw` | Processes per-trait IVW results |
| `make_ivw_parquet` | Combines all results into a single Parquet and cleans intermediate files |

Final output: `{ivw_output_dir}/.make_parquet.done`

### MR-link-2 (`all_mrlink2`)

| Rule | Description |
|------|-------------|
| `preprocess_mrlink2_outcome` | Converts GWAS `.ma` files to Parquet per trait (4 CPUs, 32 GB) |
| `preprocess_mrlink2_exposure` | Prepares eQTL exposure data per cell type (4 CPUs, 32 GB) |
| `preprocess_mrlink2_exposure_traits` | Filters gene lists per trait × cell type (4 CPUs, 32 GB) |
| `run_mrlink2` | Runs `mr_link_2_standalone.py` per gene × chromosome (8 CPUs, 32 GB, express queue) |
| `process_mrlink2` | Aggregates per-trait MR-link-2 results |
| `combine_mrlink2_results` | Merges all results into a single Parquet and cleans intermediate files |

Final output: `{mrlink2_output_dir}/results/tenk10k_phase1_sensitivity.mrlink2.parquet.gz`


## Known limitations

- **MR-link-2 is not vendored.** `rules/mrlink2.smk` invokes
  `mr_link_2_standalone.py` from a path inside the container's `/workspace` bind
  mount. That script is not in this repository, is not in the published
  container, and the scratch directory it was run from no longer exists, so the
  `all_mrlink2` target cannot currently be run by anyone. Obtain it from
  [adriaan-vd-graaf/mrlink2](https://github.com/adriaan-vd-graaf/mrlink2) and pin
  the commit you use.
- **The container has no build recipe.** `mr.sif` is referenced by every rule but
  no `.def` or `Dockerfile` exists in either repository, and the Docker Hub
  reference is untagged. Until a pinned digest or recipe is added, the software
  versions used are not recoverable.
- **`ivw_output_dir` is not fully honoured.** `rules/snakescripts/ivw_mr/process_ivw_corr.R`
  and `make_ivw_corr_parquet.R` hardcode `/workspace/ivw-ld`, so changing that
  config key moves the sentinel files but not the data the R scripts read and
  write. Change both.
- **The input has no producer.** `key_files.smr_results` points at
  `tenk10k_phase1_sensitivity.msmr.parquet.gz`, the relaxed-threshold SMR run.
  No rule in this repository produces it.
- **IVW-LD, by contrast, is complete.** `rules/snakescripts/ivw_mr/run_ivw_corr.R`
  implements it in full and takes every path as an argument.
