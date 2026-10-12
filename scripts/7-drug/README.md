# Single-cell MR supports drug target discovery and clinical translation

These scripts cross-reference the single-cell MR results with target–indication
(T–I) pairs from the Open Targets Platform (release 26.03). They estimate
enrichment of MR-supported pairs across clinical development stages and compare
Open Targets association scores between pairs with and without MR support. They
also test whether the MR direction of effect agrees with the drug mechanism of
action, and project the power of that concordance test at larger sample sizes.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `otp_combined.R` | Fig 5a–d; Extended Data Fig 8; Supp Tables 11–14 | `make fig-m5`, `make fig-e8`, `make tab-s11`, `make tab-s12`, `make tab-s13`, `make tab-s14` | shared results via `scripts/0-preprocess/preprocess_results.R`; `metadata/otp.xlsx`; `resources/misc/chembl_action_type.tsv`; `results/otp/26.03/{all_disease,all_drug_mechanism_chembl,all_drug_evidence_chembl,all_drug_molecule_chembl}.tsv`; `results/otp/26.03/otp_assoc_{overall,datasource}.gz.parquet` | `figures/fig-m5.{pdf,png}`, `figures/fig-e8.{pdf,png}`, `tables/tab-s11.tsv` … `tables/tab-s14.tsv`, `results/otp/26.03/otp_concordance_by_phase.tsv` |
| `otp_power_projection.R` | Supp Fig 16; Supp Table 15 | `make fig-s16`, `make tab-s15` | `results/otp/26.03/otp_concordance_by_phase.tsv` (from `otp_combined.R`) | `figures/fig-s16.{pdf,png}`, `tables/tab-s15.tsv` |

The Open Targets inputs under `results/otp/26.03/` come from the pipeline
(`workflow/rules/otp.smk`). All `write_table()` calls rename columns from
`metadata/table_column_names.xlsx`.

## Run order

1. `otp_combined.R` writes `results/otp/26.03/otp_concordance_by_phase.tsv`, the
   per-phase concordance between MR direction and mechanism of action.
2. `otp_power_projection.R` reads that file.

## Restricted inputs

None. Both scripts use summary-level MR results and public Open Targets data.
