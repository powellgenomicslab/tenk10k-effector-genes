# Single-cell MR supports drug target discovery and clinical translation

These scripts cross-reference the single-cell MR results with target–indication
(T–I) pairs from the Open Targets Platform (release 26.03). They estimate
enrichment of MR-supported pairs across clinical development stages and compare
Open Targets association scores between pairs with and without MR support. They
also test whether the MR direction of effect agrees with the drug mechanism of
action, and project the power of that concordance test at larger sample sizes.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `otp_combined.R` | Fig 5a–d; Supp Tables 11, 12, 13 (drugs per target–indication pair) and 15 (OT score by evidence source) | `make fig-m5`, `make tab-s11`, `make tab-s12`, `make tab-s13`, `make tab-s15` | shared results via `scripts/0-preprocess/preprocess_results.R`; `metadata/otp.xlsx`; `resources/misc/chembl_action_type.tsv`; `results/otp/26.03/{all_disease,all_drug_mechanism_chembl,all_drug_evidence_chembl,all_drug_molecule_chembl}.tsv`; `results/otp/26.03/otp_assoc_{overall,datasource}.gz.parquet` | `figures/fig-m5.{pdf,png}`, `tables/tab-s11.tsv`, `tab-s12.tsv`, `tab-s13.tsv`, `tab-s15.tsv`, `results/otp/26.03/otp_concordance_by_phase.tsv` |
| `otp_ti_mr_all.R` | Extended Data Fig 8 (association scores of MR-supported target–indication pairs by disease) | `make fig-e8` | shared results via `preprocess_results.R`; `metadata/otp.xlsx`; `results/otp/26.03/{all_disease,all_drug_evidence_chembl}.tsv`, `otp_assoc_overall.gz.parquet` | `figures/fig-e8.{pdf,png}` |
| `otp_power_projection.R` | Supp Fig 16; Supp Table 14 (MoA–MR concordance and power) | `make fig-s16`, `make tab-s14` | `results/otp/26.03/otp_concordance_by_phase.tsv` (from `otp_combined.R`) | `figures/fig-s16.{pdf,png}`, `tables/tab-s14.tsv` |

The Open Targets inputs under `results/otp/26.03/` come from the pipeline
(`workflow/rules/otp.smk`). All `write_table()` calls rename columns from
`metadata/table_column_names.xlsx`.

## Run order

1. `otp_combined.R` writes `results/otp/26.03/otp_concordance_by_phase.tsv`, the
   per-phase concordance between MR direction and mechanism of action.
2. `otp_power_projection.R` reads that file.

## Restricted inputs

None. Both scripts use summary-level MR results and public Open Targets data.
