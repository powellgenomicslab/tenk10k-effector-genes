# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is the code accompanying the TenK10K Phase 1 causal inference manuscript (medRxiv 2025.08.28.25334614). It holds Snakemake pipelines and R/Python analysis scripts only — no data, and no results.

The single most important thing to understand is that **the pipelines do not run in place.** `workflow/` and `sensitivity/` execute from the working repository at `/g/data/fy54/analysis/tenk10k-causal`, which is a separate git repo (`powellgenomicslab/tenk10k-causal`) containing the `resources/`, `results/`, `logs/`, `containers/` and `config/path/` directories that every rule and script expects. Three independent facts establish this:

- Rules in `workflow/rules/*.smk` name their inputs as `workflow/rules/snakescripts/...`, so the working directory is the *parent* of `workflow/`, not `workflow/` itself.
- Rule outputs land in `resources/`, `results/` and `logs/`, which are siblings of `workflow/` and absent from this checkout.
- `get_file_map()` in `workflow/rules/init.smk` opens `config/path/{host}.yaml` at parse time. That file exists only in `tenk10k-causal` (as `config/path/nci.yaml` and `config/path/brenner.yaml`), so `workflow/snakefile` cannot even be parsed from here.

Consequently, `workflow/README.md`'s instruction to `cd workflow/ && snakemake --profile profiles/default` contradicts the paths written into the rules; treat the documented commands as a record of intent rather than something you can execute. The invocations in `sensitivity/README.md` are likewise unverified — they pass `--profile` a YAML *file* rather than a directory, and that file mixes `configfile` keys (`key_files`, `container_paths`) with profile keys.

The figure and table scripts in `scripts/` *can* run here: `scripts/Snakefile` (front door: the root `Makefile`) declares every input, and `make link-data DATA_DIR=...` or `make data` supplies them. See "Figure and table harness" below.

Changes flow **from** `tenk10k-causal` **to** here, not the other way. This checkout is a curated subset, with the working repo's flat `scripts/` reorganised into one folder per Results subsection of the manuscript. When porting, apply the same conversions as the `v4` section-sync commits: `write_gs` → `write_table` with the Supp Table number, `resources/metadata` → `metadata/`, no `/g/data` paths or `setwd()`, read `v5`, and save figures to the standard names below. Editing a script here changes the published record, not what was actually run.

## Layout and duplication

`README.md` covers the directory structure and is current. The `scripts/` tree follows the manuscript's eight Results subsections: `0-preprocess` (shared loader and bundle helpers), `1-overview`, `2-trait`, `3-pleiotropy`, `4-comparison`, `5-pathways`, `6-polygenic`, `7-drug`, `8-crohns`, and `util`. Each folder's README lists its scripts, manuscript items and run order. The older layouts (`2-mr`, `3-comparison`, `4-polygenic`, `5-drug`, `6-crohns`, and before that `3-polygenic`, `4-drug`, `5-crohns`) are stale wherever they appear.

There is one entry point per pipeline. `sensitivity/` used to be duplicated verbatim at the repository root (`Snakefile`, `rules/`, `config/gadi.yaml`, `scripts/preparePhenotypeList.R`), which let you run `snakemake` from the root without `--snakefile`. That duplicate has been removed: the root `Snakefile` pulled in the root `rules/` but hardcoded `configfile: "sensitivity/config/gadi.yaml"`, so it ran duplicated code against the other copy's configuration, and the root `config/gadi.yaml` was never read at all. Do not reintroduce it. `config/` now holds only `config/path/{nci,brenner}.yaml`, the parse-time file maps for the main workflow.

## Parse-time configuration

Three files are read while the workflow is being *parsed*, so all of them must
exist before any rule can be evaluated. All three are now tracked:

- `config/path/nci.yaml` and `config/path/brenner.yaml` — file maps read by
  `get_file_map()` in `workflow/rules/init.smk`. They ship as templates with
  placeholder `source:` paths; a user replaces those with their own.
- `workflow/config/mvcoloc.yaml` — `configfile:` at `workflow/rules/coloc.smk:6`.
- `resources/misc/target_phenotypes.txt` — read at module scope by `_get_traits()`
  in `coloc.smk`.

`sensitivity/config/local.yaml` is the equivalent for the sensitivity pipeline
off-cluster. Verify any change with a parse check, which should exit 0:

```bash
snakemake --snakefile workflow/snakefile --list-rules
snakemake --snakefile sensitivity/Snakefile --list-rules
```

## Broken references

Every live `source()` target in `scripts/` resolves. Two `source()` calls in
`scripts/2-trait/` name files that do not exist (`scripts/mr_overview/{rrho_function,temp}.R`),
but both are inside commented-out lines. `scripts/util/helper.R` is deliberately absent — it hardcoded a private
Google Sheets ID and pulled in an OAuth dependency. Its `write_gs()` callers were
ported to `write_table()` in `scripts/util/write_table.R`, whose first three
positional arguments match. Do not reintroduce `helper.R`.

## The analysis-script contract

`scripts/0-preprocess/preprocess_results.R` is the hinge of the whole `scripts/` tree. Nearly every downstream script begins by sourcing it, and it defines the shared vocabulary:

- It reads from the tracked `metadata/` directory: `cell.tsv` (cell types, colours, ordering), `trait.tsv` (trait definitions, sample sizes, `include` flags), `trait_category.tsv` (category ordering) and `gencode.v44.gene_type.tsv.gz` (gene annotation). The older `resources/metadata/*` and `resources/misc/cell_map.tsv` paths have been migrated away; do not reintroduce them. `metadata/cell.tsv` is authoritative — the Crohn's DEG script derives the same table but must not write back over it.
- The threshold is everywhere local FDR from the `qvalue` package at `lfdr < 0.05`, but the scope differs by dataset and this matters: eQTLGen mSMR is computed **per phenotype** (`.by = "phenotype"`), MAGMA is computed **globally** across the filtered table, and TenK10K's own lfdr arrives **pre-computed** from the pipeline in `results/preprocessed/tenk10k_phase1.v5.parquet.gz` and is not recalculated here. Colocalisation is `PP.H4.abf >= 0.8`.
- Two version numbers are in play and must not be conflated: the **code** version (this branch, `v4`) and the **summary-statistics** version (`v5`). All scripts read `v5`; its changelog lives in `tenk10k-causal` (`CHANGELOG.v5.md`). The coloc vintage that goes with `v5` is `coloc.v3`, which the pipeline cannot currently regenerate — see the note at `rule concat_coloc_all`.
- The gene universe is the intersection of genes tested by both MAGMA and TenK10K MR; eQTLGen and MAGMA overlap flags are joined onto the TenK10K table with `data.table` update-joins.
- Cell type ordering comes from the `major_cell_type_order` vector defined near the top of the file. Reuse it so figures stay consistent.

Style within these scripts is a deliberate mix of `data.table` for the heavy joins and `dplyr`/`tidyverse` for reshaping, with `arrow::read_parquet` for all pipeline outputs. Match whichever idiom dominates the file you are editing.

`scripts/util/write_table.R` provides `write_table()`, the standard way to emit a manuscript supplementary table. It renames columns from `metadata/table_column_names.xlsx` (a sheet per table) and writes Supplementary Table N to `tables/tab-sN.tsv`. Every supplementary table, including Crohn's Tables 16–21, goes through it; add a sheet when adding a table.

Figures are written with `ggsave`/`ragg::agg_png` (PNG) and `cairo_pdf` (PDF) to standard names that are not committed: Main Figure N → `figures/fig-mN.{pdf,png}`, Extended Data → `figures/fig-eN.*`, Supplementary → `figures/fig-sN.*`. Panels that another script assembles go to `figures/panels/<folder>/`. Figures larger than the journal's maximum print size are shrunk after saving with `fit_journal_pdf()` (`scripts/util/fit_journal.R`, which needs `typst` and `poppler`), as in `tenk10k-causal`. Only `figures/biorender/` is tracked. Intermediate R objects (RDS/RData) stay under `results/` or `resources/crohns_case_study/`.

## Figure and table harness

`make fig-m3`, `make fig-e3`, `make fig-s3` and `make tab-s3` build single items; `make list` shows them all. The Makefile wraps `scripts/Snakefile`, which has one producer rule per script and one item rule per manuscript item. Declare every external input with `ext()` (add `review=` for cell-level or third-party files) so that `make manifest` keeps `config/figure_inputs.txt`, the bundle allow-list, complete. Steps that read raw single-cell data live under `if RESTRICTED:` and run only with `--config restricted=true`. Steps that read individual-level genotypes belong in `workflow/` (the GPX1 and JAZF1 locus objects are built by rules in `workflow/rules/misc.smk`); the harness takes their summary-level outputs (LD r² only, never genotypes) as `ext()` inputs. Software comes from `workflow/envs/figures.yaml`, pinned per platform in `figures.<platform>.pin.txt` (`make lock`); `make dockerfile` regenerates the Dockerfile with the pins applied. Check changes with `make dry-run` and `make list`. The Zenodo bundle is built with `make bundle` (after `make link-data`): `make inventory` first writes `INVENTORY.tsv`, joining each manifest file with its size, SHA-256, the display items that read it (traced through the Snakefile) and the description, source and licence from `config/figure_inputs_annotation.tsv`, whose first matching pattern wins. Add a row there when adding an input; `make inventory` fails on any file without one. The archive root also carries `scripts/0-preprocess/BUNDLE_README.md` as `README.md`.

## Snakemake specifics

Two pipelines, both driven from the working repository root:

- `workflow/snakefile` is the main pipeline: GWAS formatting and liftover, BESD/SMR preparation, SMR and SMR-multi, MAGMA, coloc and multi-variant coloc, scDRS, genetic correlation via LDAK, gene-set enrichment (Enrichr/g:Profiler/STRING), and Open Targets Platform queries. Rule files live in `workflow/rules/*.smk`, with the code they call in `workflow/rules/snakescripts/`.
- `sensitivity/Snakefile` runs the sensitivity MR methods, IVW-LD and MR-link-2, inside the `mr.sif` Singularity container. It reads the main pipeline's SMR-multi Parquet output, so the main pipeline must have completed first. Targets are `all_ivw` and `all_mrlink2`; a checkpoint, `prepare_phenotype_list`, resolves the trait list dynamically.

Conventions to follow when adding rules:

- A rule that dispatches to a per-study script does so through `params.script` with a wildcard in the path, e.g. `params: script = "snakescripts/prep_smr_input/{study}.R"` then `script: "{params.script}"`. Add a new study by dropping in a file named after it, not by editing the rule.
- Conda environments are defined in `workflow/envs/{renv,scverse,pydata,sc-renv}.yaml` and referenced from rules by relative path (`conda: "../envs/renv.yaml"`, resolved relative to the `.smk` file). Versions were captured from the environment that produced the published results (R 4.4.1). `sc-renv` additionally needs four GitHub-only packages installed by hand; its header says which. The top-level `environment.yml` describes an unrelated environment and is not what the pipeline uses.
- `workflow/profiles/default/` contains only `config.v8+.yaml`, with no `config.yaml` fallback, so the main pipeline requires Snakemake 8 or newer despite `environment.yml` asking for `>=7`. It submits to PBS on NCI Gadi via `cluster-generic`, project `fy54`, and needs `-l storage` to cover `gdata/ei56`, `gdata/fy54` and the matching scratch filesystems.
- Sensitivity rules run under Singularity with fixed bind mounts (`/data`, `/genotypes`, `/workspace`, `/eqtl`); container-internal paths are configured separately under `container_paths` in `sensitivity/config/gadi.yaml`, distinct from the host paths under `key_files`. When adding a rule there, both must be set. Two shell blocks in `sensitivity/rules/mrlink2.smk` are raw strings (`r"""`) because the shell must receive `\t` and `\;` verbatim; keep them raw.
- A good deal of commented-out rule code is retained in the `.smk` files as a record of superseded approaches. Leave it be unless asked.

Snakemake is not on `PATH` in a bare login shell on Gadi; it comes from a module or conda environment that must be loaded first.

## Housekeeping

There are no tests and no linter configured — do not invent commands for either.

`.gitignore` blanket-ignores `*.tsv`, `*.csv`, `*.xlsx`, `*.parquet`, `*.gz` and similar at any depth, as well as `results/`, `data/` and `logs/`. The tracked files under `metadata/` were therefore force-added, and any new data or metadata file will need `git add -f`. Check what is actually tracked with `git ls-files` rather than assuming.

The `TenK10K_causal_inference_vs_IBDverse` submodule is not checked out; run `git submodule update --init` if you need it.

Work happens on the `v4` branch; `main` is the default for pull requests.
