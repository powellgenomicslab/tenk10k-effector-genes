
# extract gwas locus

rule preprocess_results:
    """"
    Preprocess results for downstream analysis
    """
    output:
        "results/preprocessed/{study}.{version}.parquet.gz"
    conda: "../envs/renv.yaml"
    resources:
        mem = "64GB",
        ncpus = 16
    script:
        "snakescripts/preprocess/{wildcards.study}.{wildcards.version}.R"

rule extract_gwas_locus:
    input:
        gwas = "resources/ma/{pheno}.ma",
        bfile = expand("resources/genotypes/{{study}}/chr{{chr}}.{ext}", ext = ["bed", "bim", "fam"])
    output:
        ld = "results/gwas_locus/{pheno}/{study}/{chr}_{start}_{end}.ld",
        gwas = "results/gwas_locus/{pheno}/{study}/{chr}_{start}_{end}.gwas.tsv"
    params:
        bfile_prefix = "resources/genotypes/{study}/chr{chr}",
        out_prefix = "results/gwas_locus/{pheno}/{study}/{chr}_{start}_{end}"
    log:
        "logs/gwas_locus/{pheno}/{study}/{chr}_{start}_{end}.log"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.gwas}) && \
        awk -v OFS='\\t' '
            NR==1 {{print "chr", "pos", $0; next}}
            {{split($1, a, ":")
              if (a[1] == "{wildcards.chr}" && a[2] >= {wildcards.start} && a[2] <= {wildcards.end})
              print a[1], a[2], $0}}' \
            {input.gwas} > {output.gwas}
        
        TOPVAR=$(awk 'NR == FNR {{a[$2];next}}
         $3 in a {{if (p == "" || $9 < p) {{p = $9; var = $3}}}}
         END {{print var}}' {input.bfile[1]} {output.gwas})
        
        echo "Top variant: $TOPVAR"
        
        plink --bfile {params.bfile_prefix} \
              --chr {wildcards.chr} \
              --from-bp {wildcards.start} \
              --to-bp {wildcards.end} \
              --silent \
              --ld-window-r2 0 \
              --ld-window 1e6 \
              --r2 \
              --ld-snp $TOPVAR \
              --out {params.out_prefix}
        """
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
        "Rscript workflow/rules/snakescripts/crohns_locus/1-generate_locus_zoom_object.R "
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
        "Rscript workflow/rules/snakescripts/crohns_locus/locus_zoom_eQTL_discordant.R > {log} 2>&1"

# Extracts of two published supplementary tables (Kong et al. 2023, Liu et al. 2023)
# used by the Crohn's case study. The originals are not redistributed; the figure
# harness reads these extracts, which are shared in the data bundle with a citation.

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
        "Rscript workflow/rules/snakescripts/crohns_external/extract_third_party.R > {log} 2>&1"
