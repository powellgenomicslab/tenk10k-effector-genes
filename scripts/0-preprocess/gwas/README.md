# GWAS summary statistics preparation

The code that harmonises GWAS summary statistics does **not** live here — it is
part of the Snakemake pipeline. This directory is kept as a signpost.

| Step | Location |
|---|---|
| Per-trait formatting into GCTA `.ma` format (25 traits) | `workflow/rules/snakescripts/format_gwas/{trait}.R`, driven by `rule format_gwas` |
| GRCh37 → GRCh38 liftover helper, sourced by 18 of those scripts | `workflow/rules/snakescripts/hg19tohg38.R` |
| FinnGen/MVP/UKBB meta-analysis extraction | `rule extract_finngen_gwas` + `rule format_finngen_gwas`, configured by `workflow/config/finngen_meta_path.yaml` |
| Splitting `.ma` files by chromosome for coloc | `rule split_ma_by_chr` |

Trait definitions, sample sizes and inclusion flags are in `metadata/trait.tsv`.

## How `resources/ma/` is assembled

Every downstream arm (SMR, MAGMA, coloc, LDAK and the sensitivity MR pipeline)
reads `resources/ma/{trait}.ma`. That directory is a curated **input**, not a
pipeline output, and no rule writes it. It is assembled from three routes:

1. `rule format_gwas` writes harmonised files to the `resources/pipeline_ma/`
   staging area; these are reviewed before being promoted into `resources/ma/`.
2. The FinnGen route writes to `resources/sumstats/finngen_gwas_extract/`.
3. Summary statistics already supplied in `.ma` format by the source consortium
   are copied in directly.

Anyone reproducing the analysis must assemble `resources/ma/` the same way. See
the provenance note at the top of `workflow/rules/format_gwas.smk`.
