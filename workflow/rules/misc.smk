
# extract gwas locus

rule preprocess_results:
    """"
    Preprocess results for downstream analysis
    """
    output:
        "results/preprocessed/{study}.{version}.parquet.gz"
    conda: "renv"
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