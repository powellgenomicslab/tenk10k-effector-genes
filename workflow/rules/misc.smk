
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
# The main results as a standalone file for the Zenodo data record: every column
# of the integrated table, with clear names and with the column descriptions of
# workflow/config/main_results_columns.tsv stored in the Parquet metadata.

rule main_results_release:
    input:
        results = "results/preprocessed/tenk10k_phase1.v5.parquet.gz"
    output:
        parquet = "results/release/tenk10k_phase1_effector_genes.v5.parquet",
        md = "results/release/tenk10k_phase1_effector_genes.v5.md"
    params:
        columns = f"{workflow.basedir}/config/main_results_columns.tsv",
        meta = {
            "title": "TenK10K Phase 1 effector genes: main results",
            "description": (
                "Single-cell eQTL Mendelian randomisation results for every gene, cell type "
                "and trait tested in TenK10K Phase 1: multi-SNP SMR with local FDR, the "
                "HEIDI, IVW-LD and MR-link-2 sensitivity analyses, coloc and SuSiE-based "
                "colocalisation, and the evidence tiers used in the manuscript."
            ),
            "version": "summary statistics v5, colocalisation v3",
            "doi": "10.5281/zenodo.23406585",
            "licence": "CC BY-SA 4.0",
            "citation": (
                "Henry A, Senabouth A, Tyebally R, et al. Single-cell genetics identifies "
                "cell type-specific causal mechanisms in complex traits and diseases. "
                "medRxiv 2025.08.28.25334614 (2025). doi:10.1101/2025.08.28.25334614"
            ),
            "code": "https://github.com/powellgenomicslab/tenk10k-effector-genes",
        }
    log:
        "logs/release/main_results.log"
    conda: "scverse"
    script: "snakescripts/release/main_results.py"
