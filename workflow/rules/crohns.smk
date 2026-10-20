# Crohn's disease case study: raw and third-party inputs -> the intermediate files
# shared in the data bundle (Fig. 6, ED Fig. 9, Supp. Figs 17 and 25, Supp. Tables
# 16-21). The figure harness (scripts/Snakefile) takes these outputs from the bundle.
# Some rules read MR results annotated for Crohn's disease
# (resources/crohns_case_study/postprocess/), which the harness's
# scripts/8-crohns/prepare_data/2-postprocess.R writes.

# Crohn's disease locus objects for the manuscript figures (ED Fig. 9, Supp. Fig. 17).
# Each combines eQTL and GWAS statistics with LD (r2 to the index variant) computed
# by PLINK from the TenK10K genotypes; only the r2 values are written out, never
# the genotypes. PLINK must be on PATH, as for extract_gwas_locus.

rule crohns_gpx1_locus_data:
    input:
        gwas = "resources/pipeline_ma/crohns.ma",
        eqtl = expand("resources/matrix_eqtl/tenk10k_phase1/{ct}/chr3_meqtl.tsv", ct = ["CD4_TCM", "cDC2"]),
        bfile = expand("resources/genotypes/tenk10k_phase1/chr3.{ext}", ext = ["bed", "bim", "fam"])
    output:
        "resources/crohns_case_study/locus_zoom_data/GPX1/tenk10k_GPX1_locus_plot_data.RDS"
    log:
        "logs/crohns_locus/GPX1.log"
    conda: "../envs/figures.yaml"
    shell:
        "Rscript workflow/rules/snakescripts/crohns/1-generate_locus_zoom_object.R "
        "GPX1 3 ENSG00000233276 > {log} 2>&1"

rule crohns_jazf1_locus_data:
    input:
        # MR results annotated for Crohn's disease, written by
        # scripts/8-crohns/prepare_data/2-postprocess.R; only the JAZF1 probe ID,
        # chromosome and top eQTL variants are read from it
        mr = "resources/crohns_case_study/postprocess/tenk_crohns_raw.RDS",
        eqtl = expand("resources/matrix_eqtl/tenk10k_phase1/{ct}/chr7_meqtl.tsv", ct = ["CD4_Naive", "cDC2"]),
        bfile = expand("resources/genotypes/tenk10k_phase1/chr7.{ext}", ext = ["bed", "bim", "fam"])
    output:
        "resources/crohns_case_study/figures/JAZF1_data.RDS"
    log:
        "logs/crohns_locus/JAZF1.log"
    conda: "../envs/figures.yaml"
    shell:
        "mkdir -p resources/crohns_case_study/figures/ld && "
        "Rscript workflow/rules/snakescripts/crohns/locus_zoom_eQTL_discordant.R > {log} 2>&1"

# Extracts of two published supplementary tables (Kong et al. 2023, Liu et al. 2023)
# used by the Crohn's case study, read by crohns_deg_preprocess and
# crohns_external_annot below. Neither the originals nor the extracts are shared.

rule crohns_third_party_extracts:
    input:
        kong = "resources/crohns_case_study/deg/Kongetal2023_supplementary/Kongetal2023_Crohns-sc_DEG.xlsx",
        features = "resources/crohns_case_study/deg/features_original.tsv",
        liu = "resources/crohns_case_study/crohns_relevant_gene_lists/sources/41588_2023_1384_MOESM4_ESM_Liuetal2023_SupTable.xlsx"
    output:
        kong = "resources/crohns_case_study/deg/kong2023_immune_deg.parquet",
        liu = "resources/crohns_case_study/crohns_relevant_gene_lists/sources/liu2023_st8_genes.tsv"
    log:
        "logs/crohns_locus/third_party_extracts.log"
    conda: "../envs/figures.yaml"
    shell:
        "Rscript workflow/rules/snakescripts/crohns/extract_third_party.R > {log} 2>&1"

# Literature evidence per gene: membership of the Crohn's disease gene lists
# (systematic review, Open Targets, Liu et al. 2023, drug targets and pathways).

rule crohns_external_annot:
    input:
        gene_annot = "resources/misc/gencode.v44.gene_type.tsv",
        lists = expand("resources/crohns_case_study/crohns_relevant_gene_lists/sources/{f}", f = [
            "Systematic Review of Crohns Disease Genes.csv",
            "OT-EFO_0000384-associated-targets-07_04_2026-v26_03.tsv",
            "liu2023_st8_genes.tsv",
            "drug_target_genes.txt",
        ])
    output:
        "resources/crohns_case_study/crohns_relevant_gene_lists/external_source_gene_annot_df.RDS"
    log:
        "logs/crohns/external_annot.log"
    conda: "../envs/figures.yaml"
    shell:
        "Rscript workflow/rules/snakescripts/crohns/1-prepare_external_source_annot_boolean_df.R > {log} 2>&1"

# Single-cell DEGs of Kong et al. 2023 by major cell type, joined with the Crohn's
# disease MR genes; and the number of DEG genes in the matched major cell types,
# which Fig. 6a reports.

rule crohns_deg_preprocess:
    input:
        cell_map = "resources/metadata/cell_map_revised.tsv",
        deg = "resources/crohns_case_study/deg/kong2023_immune_deg.parquet"
    output:
        "resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS",
        "resources/crohns_case_study/deg/common_gut_cell_types_with_tenk_revision.txt"
    log:
        "logs/crohns/deg_preprocess.log"
    conda: "../envs/figures.yaml"
    shell:
        "Rscript workflow/rules/snakescripts/crohns/1-prepare_crohns_deg_pre-processed.R > {log} 2>&1"

rule crohns_mr_deg:
    input:
        "resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS",
        "resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS"
    output:
        "resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS",
        expand("resources/crohns_case_study/deg/mr_deg_{x}.txt", x = ["con", "discon", "mixed"])
    log:
        "logs/crohns/mr_deg.log"
    conda: "../envs/figures.yaml"
    shell:
        "Rscript workflow/rules/snakescripts/crohns/2-mr_deg_comparison_inner_join.R > {log} 2>&1"

rule crohns_deg_count:
    input:
        "resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS",
        "resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS"
    output:
        "resources/crohns_case_study/deg/deg_genes_matched_major_ct.tsv"
    log:
        "logs/crohns/deg_count.log"
    conda: "../envs/figures.yaml"
    shell:
        "Rscript workflow/rules/snakescripts/crohns/5-count_matched_deg_genes.R > {log} 2>&1"

# Supp. Fig. 25: harmonised major cell types on the UMAP of the Kong et al. 2023
# colon immune cells. It needs their single-cell data, which is not shared, so this
# figure is drawn here and cannot be rebuilt from the data bundle.

rule crohns_h5ad_major_ct:
    input:
        "resources/crohns_case_study/deg/cd_colon_immune.h5ad"
    output:
        "resources/crohns_case_study/deg/cd_colon_immune_major_cell_types.h5ad"
    log:
        "logs/crohns/h5ad_major_ct.log"
    conda: "scverse"
    shell:
        "python3 workflow/rules/snakescripts/crohns/3-prepare_data_for_MR_comparison.py > {log} 2>&1"

rule crohns_umap:
    input:
        "resources/crohns_case_study/deg/cd_colon_immune_major_cell_types.h5ad"
    output:
        expand("figures/panels/8-crohns/umap_{x}.{e}", x = ["original_annotation", "major_cell_types"], e = ["png", "pdf"])
    log:
        "logs/crohns/umap.log"
    conda: "scverse"
    shell:
        "python3 workflow/rules/snakescripts/crohns/4-plot_UMAP_cell_annotation.py > {log} 2>&1"
