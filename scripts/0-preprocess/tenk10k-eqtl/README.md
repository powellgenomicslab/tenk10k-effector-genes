# TenK10K Phase 1 sc-eQTL preparation

The sc-eQTL mapping itself is a separate study (Cuomo et al., medRxiv
2025.03.20.25324352) and is **not** reproduced here. Nor does this directory hold
code: the steps that bring SAIGE-QTL output into this pipeline are Snakemake
rules. This directory is kept as a signpost.

| Step | Location |
|---|---|
| Stage pre-built BESD/ESI files, build per-chromosome EPI files | `workflow/rules/snakescripts/prep_besd_chr/tenk10k_phase1.sh` and `prep_epi_chr.tenk10k_phase1.R`, driven by `rule prep_besd_chr` |
| SMR probe lists and per-cell-type p-value thresholds | `workflow/rules/snakescripts/prep_smr_input/tenk10k_phase1.R` |
| Per-chromosome eQTL input for coloc (`.fst`) | `workflow/rules/snakescripts/prep_coloc_input/tenk10k_phase1.R` |
| eQTL LD reference | `workflow/rules/snakescripts/make_ld/eqtl/tenk10k_phase1.R` |

## What you must supply

- **Pre-built BESD files.** `prep_besd_chr` stages them; it does not build them.
  The SAIGE-QTL → BESD conversion lives in the external `TenK10K_SMR`
  repository. Set `TENK10K_BESD_DIR` and `GENCODE_GTF` to point at your copies.
- **`resources/brenner/tenk10k_phase1/common_eqtl.tsv`** — the common eQTL /
  eGene table from the sc-eQTL study. `config/path/brenner.yaml` documents where
  it comes from and symlinks it into place.
- **`resources/saige_eqtl/tenk10k_phase1/{cell_type}/common_raw.tsv`** — raw
  per-cell-type association statistics, required by the coloc input step.
