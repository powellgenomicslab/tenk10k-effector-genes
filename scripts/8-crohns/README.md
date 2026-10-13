# Single-cell genetics coupled with triangulation of evidence prioritises effector genes for Crohn’s disease

These scripts carry out the Crohn's disease case study. They take the Crohn's MR associations that pass the sensitivity and colocalisation filters and annotate them with several sources:

- literature and drug-target gene lists;
- bulk eQTLGen MR and MAGMA;
- differential expression in intestinal biopsies (Kong et al. 2023).

They then compare the results with IBDverse. Together the scripts produce Figure 6a–c, Extended Data Figure 9, Supplementary Figures 17, 25 and 30, and Supplementary Tables 16–21. All scripts run from the repository root.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `prepare_data/1-prepare_external_source_annot_boolean_df.R` | Supp Tab 21 (input) | – | `resources/crohns_case_study/crohns_relevant_gene_lists/sources/*`, `metadata/gencode.v44.gene_type.tsv.gz` | `resources/crohns_case_study/crohns_relevant_gene_lists/external_source_gene_annot_df.RDS` |
| `prepare_data/2-postprocess.R` | all (input) | – | `scripts/0-preprocess/preprocess_results.R` inputs, coloc v3 and mvcoloc parquet, `external_source_gene_annot_df.RDS` | `resources/crohns_case_study/postprocess/{tenk_crohns_raw,tenk_crohns_sig,tenk_alltraits_raw,tenk_alltraits_sig,eqtlgen_alltraits_raw,df_coloc,df_mvcoloc}.RDS` |
| `deg/1-prepare_crohns_deg_pre-processed.R` | Fig 6a–c, Supp Tab 18 (input) | – | Kong et al. 2023 DEG supplementary xlsx, `resources/crohns_case_study/deg/features_original.tsv`, `metadata/cell.tsv` | `resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS`, `common_gut_cell_types_with_tenk_revision.txt` |
| `deg/2-mr_deg_comparison_inner_join.R` | Fig 6a–c, Supp Tab 16, 18 (input) | – | `crohns_deg_pre-processed_revision.RDS`, `tenk_crohns_sig.RDS` | `resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS`, `mr_deg_{con,discon,mixed}.txt` |
| `deg/3-prepare_data_for_MR_comparison.py` | Supp Fig 30 (input) | – | `resources/crohns_case_study/deg/cd_colon_immune.h5ad` | `resources/crohns_case_study/deg/cd_colon_immune_major_cell_types.h5ad` |
| `deg/4-plot_UMAP_cell_annotation.py` | Supp Fig 30 | `make fig-s30` | `cd_colon_immune_major_cell_types.h5ad` | `figures/panels/8-crohns/umap_{original_annotation,major_cell_types}.{png,pdf}` |
| `plot_figures/1-gene_annotation_barplot.R` | Fig 6a | `make fig-m6` | `tenk_crohns_sig.RDS`, `mr_max_evidence_and_deg_innerjoin_by_majct.RDS`, `external_source_gene_annot_df.RDS`, `crohns_deg_pre-processed_revision.RDS`, `eqtlgen_alltraits_raw.RDS` | `figures/panels/8-crohns/fig6a_ext_source_annotation_barplot.{png,pdf}` |
| `plot_figures/3-annotated_heatmap.R` | Fig 6b | `make fig-m6` | `tenk_crohns_raw.RDS`, `mr_max_evidence_and_deg_innerjoin_by_majct.RDS` | `figures/panels/8-crohns/fig6b_heatmap_canonical_deg_mrsenscoloc.{png,pdf}`, `resources/crohns_case_study/figures/heatmap_objects.RData` |
| `plot_figures/4-mr_deg_scatterplot.R` | Fig 6c | `make fig-m6` | `mr_max_evidence_and_deg_innerjoin_by_majct.RDS` | `figures/panels/8-crohns/fig6c_mr_deg_annot.{png,pdf}`, `resources/crohns_case_study/figures/mr_deg_scatterplot_objects.RData` |
| `plot_figures/1-generate_locus_zoom_object.R` | ED Fig 9 (input) | – | `resources/pipeline_ma/crohns.ma`, `resources/matrix_eqtl/tenk10k_phase1/{CD4_TCM,cDC2}/chr3_meqtl.tsv`, `resources/genotypes/tenk10k_phase1/chr3.*` | `resources/crohns_case_study/locus_zoom_data/GPX1/tenk10k_GPX1_locus_plot_data.RDS` |
| `plot_figures/2-plot_locus_zoom_fig.R` | ED Fig 9: a) GPX1 locus zoom; b, c) eQTL vs GWAS effects of the MR instruments in CD4 TCM and cDC2 | `make fig-e9` | `tenk10k_GPX1_locus_plot_data.RDS`, `results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.snps4msmr.parquet.gz`, `results/smr_locus/tenk10k_phase1/{CD4_TCM,cDC2}/crohns/ENSG00000233276.{eqtl,gwas}.tsv` | `figures/fig-e9.{png,pdf}` |
| `plot_figures/locus_zoom_eQTL_discordant.R` | Supp Fig 17 (input) | – | `tenk_crohns_raw.RDS`, `resources/matrix_eqtl/tenk10k_phase1/{CD4_Naive,cDC2}/chr7_meqtl.tsv`, `resources/genotypes/tenk10k_phase1/chr7.*` | `resources/crohns_case_study/figures/JAZF1_data.RDS` |
| `plot_figures/JAZF1_scatterplot.R` | Supp Fig 17 | `make fig-s17` | `JAZF1_data.RDS` | `figures/fig-s17.{png,pdf}` |
| `plot_figures/ibdverse_novelty_by_celltype.R` | Supp Fig 25 | `make fig-s25` | `results/pi1_mr_novelty_results/tenk_not_novel_res_relative_celltype.tsv` (from the IBDverse comparison) | `figures/fig-s25.{png,pdf}` |
| `supplementary/1-prepare_supplementary_table_overall_results.R` | – (input to Supp Tab script) | – | `tenk_crohns_sig.RDS` | `resources/crohns_case_study/supplementary/crohns_summary_gene_count_df.RData` |
| `supplementary/2-*.R` | Supp Tab 16–18, 21 | sourced by `3-` | as above | data frames in the session |
| `supplementary/3-prepare_supplementary_tables_all.R` | Supp Tab 16–21 | `make tab-s16` … `make tab-s21` | outputs above, `resources/crohns_case_study/revision/{IBDverse_pi1_tables.RData,OR_IBDverse_coloc_genes_in_CD_MR_genes.RDS}`, `metadata/table_column_names.xlsx` | `tables/tab-s16.tsv` … `tables/tab-s21.tsv` |

## Run order

1. `prepare_data/1-…`, then `prepare_data/2-postprocess.R`.
2. `deg/1-…`, then `deg/2-…`.
   - Separately, `deg/3-…` then `deg/4-…` for Supp Fig 30.
3. `plot_figures/1-gene_annotation_barplot.R`, `3-annotated_heatmap.R` and `4-mr_deg_scatterplot.R`, for Fig 6a–c.
4. `plot_figures/1-generate_locus_zoom_object.R GPX1 3 ENSG00000233276`, then `2-plot_locus_zoom_fig.R`, for ED Fig 9.
5. `plot_figures/locus_zoom_eQTL_discordant.R`, then `JAZF1_scatterplot.R`, for Supp Fig 17.
6. `supplementary/1-…`, then `supplementary/3-…`. The second script sources the `2-*` scripts.

Figure 6 is assembled from the panels in `figures/panels/8-crohns/`, and no script does the assembly. Panel 6d is not produced here. It comes from the IBDverse comparison (`TenK10K_causal_inference_vs_IBDverse/scripts/enrichment.r`), as does Supplementary Figure 24 (π1 replication, `pi1_estimate.r`).

## Restricted inputs

- `deg/3-…` and `deg/4-…` read the Crohn's colon single-cell h5ad (Kong et al. 2023). That is raw single-cell data, and it is not shared.
- `plot_figures/1-generate_locus_zoom_object.R` and `plot_figures/locus_zoom_eQTL_discordant.R` compute LD with PLINK from individual-level TenK10K genotypes. Neither can be rerun from the public release.
- The downstream figure and table scripts read only the summary-level intermediates listed above.
