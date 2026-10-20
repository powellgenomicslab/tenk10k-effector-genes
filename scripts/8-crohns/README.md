# Single-cell genetics coupled with triangulation of evidence prioritises effector genes for Crohn’s disease

These scripts carry out the Crohn's disease case study. They take the Crohn's MR associations that pass the sensitivity and colocalisation filters and annotate them with several sources:

- literature and drug-target gene lists;
- bulk eQTLGen MR and MAGMA;
- differential expression in intestinal biopsies (Kong et al. 2023).

They then compare the results with IBDverse. The scripts in this folder turn the intermediate files of the data bundle into Figure 6, Extended Data Figures 9–10, Supplementary Figure 17 and Supplementary Tables 16–21. The steps that build those intermediate files from raw or unshared data are in `workflow/rules/crohns.smk` (see below), as is Supplementary Figure 25. All scripts run from the repository root.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `prepare_data/2-postprocess.R` | all (input) | – | `scripts/0-preprocess/preprocess_results.R` inputs, coloc v3 and mvcoloc parquet, `external_source_gene_annot_df.RDS` | `resources/crohns_case_study/postprocess/{tenk_crohns_raw,tenk_crohns_sig,tenk_alltraits_raw,tenk_alltraits_sig,eqtlgen_alltraits_raw,df_coloc,df_mvcoloc}.RDS` |
| `plot_figures/1-gene_annotation_barplot.R` | Fig 6a | `make fig-m6` | `tenk_crohns_sig.RDS`, `mr_max_evidence_and_deg_innerjoin_by_majct.RDS`, `external_source_gene_annot_df.RDS`, `deg_genes_matched_major_ct.tsv`, `eqtlgen_alltraits_raw.RDS` | `figures/panels/8-crohns/fig6a_ext_source_annotation_barplot.{png,pdf}` |
| `plot_figures/3-annotated_heatmap.R` | Fig 6b | `make fig-m6` | `tenk_crohns_raw.RDS`, `mr_max_evidence_and_deg_innerjoin_by_majct.RDS` | `figures/panels/8-crohns/fig6b_heatmap_canonical_deg_mrsenscoloc.{png,pdf}`, `resources/crohns_case_study/figures/heatmap_objects.RData` |
| `plot_figures/4-mr_deg_scatterplot.R` | Fig 6c | `make fig-m6` | `mr_max_evidence_and_deg_innerjoin_by_majct.RDS` | `figures/panels/8-crohns/fig6c_mr_deg_annot.{png,pdf}`, `resources/crohns_case_study/figures/mr_deg_scatterplot_objects.RData` |
| `plot_figures/2-plot_locus_zoom_fig.R` | ED Fig 9: a) GPX1 locus zoom; b, c) eQTL vs GWAS effects of the MR instruments in CD4 TCM and cDC2 | `make fig-e9` | `tenk10k_GPX1_locus_plot_data.RDS`, `results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.snps4msmr.parquet.gz`, `results/smr_locus/tenk10k_phase1/{CD4_TCM,cDC2}/crohns/ENSG00000233276.{eqtl,gwas}.tsv` | `figures/fig-e9.{png,pdf}` |
| `plot_figures/JAZF1_scatterplot.R` | Supp Fig 17 | `make fig-s17` | `JAZF1_data.RDS` | `figures/fig-s17.{png,pdf}` |
| `plot_figures/ibdverse_novelty_by_celltype.R` | Extended Data Fig 10 | `make fig-e10` | `results/pi1_mr_novelty_results/tenk_not_novel_res_relative_celltype.tsv` (from the IBDverse comparison) | `figures/fig-e10.{png,pdf}` |
| `supplementary/1-prepare_supplementary_table_overall_results.R` | – (input to Supp Tab script) | – | `tenk_crohns_sig.RDS` | `resources/crohns_case_study/supplementary/crohns_summary_gene_count_df.RData` |
| `supplementary/2-*.R` | Supp Tab 16–18, 21 | sourced by `3-` | as above | data frames in the session |
| `supplementary/3-prepare_supplementary_tables_all.R` | Supp Tab 16–21 | `make tab-s16` … `make tab-s21` | outputs above, `results/pi1_mr_novelty_results/{pi1_combined.tsv,pi1_all.txt,tenk_not_novel_res_relative.tsv,prop_blood_ibdverse.csv}` (from the IBDverse comparison), `metadata/cell.tsv`, `metadata/table_column_names.xlsx` | `tables/tab-s16.tsv` … `tables/tab-s21.tsv` |


## Built by the analysis pipeline (`workflow/rules/crohns.smk`)

These steps read raw or unshared data, so they run in the working repository with the
analysis pipeline, and their outputs are shared in the data bundle. Scripts are in
`workflow/rules/snakescripts/crohns/`.

| Rule | Script | Reads | Writes (shared in the bundle unless noted) |
|---|---|---|---|
| `crohns_third_party_extracts` | `extract_third_party.R` | Kong et al. 2023 DEG table and feature list; Liu et al. 2023 Supplementary Table 8 | `kong2023_immune_deg.parquet`, `liu2023_st8_genes.tsv` (not shared) |
| `crohns_external_annot` | `1-prepare_external_source_annot_boolean_df.R` | the literature gene lists, GENCODE v44 | `crohns_relevant_gene_lists/external_source_gene_annot_df.RDS` |
| `crohns_deg_preprocess` | `1-prepare_crohns_deg_pre-processed.R` | `kong2023_immune_deg.parquet`, the cell map | `deg/crohns_deg_pre-processed_revision.RDS` (not shared) |
| `crohns_mr_deg` | `2-mr_deg_comparison_inner_join.R` | the pre-processed DEGs, `tenk_crohns_sig.RDS` | `deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS` |
| `crohns_deg_count` | `5-count_matched_deg_genes.R` | the pre-processed DEGs, `tenk_crohns_sig.RDS` | `deg/deg_genes_matched_major_ct.tsv` (the count in Fig 6a) |
| `crohns_gpx1_locus_data` | `1-generate_locus_zoom_object.R` | Crohn's GWAS, CD4 TCM and cDC2 eQTL, TenK10K genotypes (chr 3) | `locus_zoom_data/GPX1/tenk10k_GPX1_locus_plot_data.RDS` (ED Fig 9) |
| `crohns_jazf1_locus_data` | `locus_zoom_eQTL_discordant.R` | `tenk_crohns_raw.RDS`, CD4 Naive and cDC2 eQTL, TenK10K genotypes (chr 7) | `figures/JAZF1_data.RDS` (Supp Fig 17) |
| `crohns_h5ad_major_ct`, `crohns_umap` | `3-prepare_data_for_MR_comparison.py`, `4-plot_UMAP_cell_annotation.py` | Kong et al. 2023 colon immune single-cell data | Supplementary Figure 25 (not rebuildable from the bundle) |

The locus objects contain eQTL and GWAS statistics and the r² of each variant with the
index variant, never genotypes. `tenk_crohns_sig.RDS` and `tenk_crohns_raw.RDS` are
written by `prepare_data/2-postprocess.R` in this folder.

## Run order

`make` (`scripts/Snakefile`) runs these in order; by hand:

1. `prepare_data/2-postprocess.R`.
2. `plot_figures/1-gene_annotation_barplot.R`, `3-annotated_heatmap.R`, `4-mr_deg_scatterplot.R` and `ibdverse_blood_enrichment.R`, for Fig 6a–d.
3. `plot_figures/2-plot_locus_zoom_fig.R`, for ED Fig 9.
4. `plot_figures/JAZF1_scatterplot.R`, for Supp Fig 17.
5. `supplementary/1-…`, then `supplementary/3-…`. The second script sources the `2-*` scripts.

Figure 6 is assembled from the panels in `figures/panels/8-crohns/` by `fig6.typ` (`make fig-m6`, compiled with Typst like Figure 1). Panel 6d is drawn by `plot_figures/ibdverse_blood_enrichment.R`, ported from the IBDverse comparison (`TenK10K_causal_inference_vs_IBDverse/scripts/enrichment.r`, branch `new_annot`), from two tables written by that analysis's `summarise.r` (`results/pi1_mr_novelty_results/{prop_blood_ibdverse.csv,tenk_not_novel_res_relative.tsv}`). Supplementary Figure 24 (π1 replication, `pi1_estimate.r`) also comes from that analysis.

## Data not shared

The literature gene lists (Garza-Hernández et al. 2022, Open Targets, Liu et al. 2023,
drug-target reviews), the single-cell DEG table of Kong et al. 2023, *Immunity* 56:444
(doi:10.1016/j.immuni.2023.01.002), and their colon single-cell data are not shared;
only the results derived from them by the pipeline above are. The figures and tables
built from those shared results are identical to those built from the sources.
