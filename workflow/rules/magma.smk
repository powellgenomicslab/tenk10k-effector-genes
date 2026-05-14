# workflow to run magma
# URL: https://cncr.nl/research/magma/


rule magma_gene_annot:
    input:
        snp_loc = "resources/genotypes/{geno_set}/chr{chr}.bim",
        gene_loc = "resources/magma/geneanno.loc"
    output:
        "results/magma/{geno_set}/chr{chr}.genes.annot"
    params:
        gene_windows_bp=0,
        prefix_out = lambda x: f"results/magma/{x.geno_set}/chr{x.chr}"
    shell:
        """
        magma --annotate window={params.gene_windows_bp} \
            --snp-loc {input.snp_loc} \
            --gene-loc {input.gene_loc} \
            --out {params.prefix_out}
        """

def get_n_trait(x):
    """
    Get the number of samples for a given pheno from metadata file
    """
    df_meta = pd.read_csv("resources/metadata/trait_metadata_n.tsv", encoding_errors = "ignore", sep="\t")
    n = df_meta.loc[df_meta['trait_id'] == x.pheno, 'n_eff'].values[0]
    return int(n)

rule magma_calc_gene_assoc:
    input:
        gwas_file = "resources/ma/{pheno}.ma",
        bfile = expand("resources/genotypes/{{geno_set}}/chr{{chr}}.{ext}", \
                        ext = ["bed", "bim", "fam"]),
        genes_annot = rules.magma_gene_annot.output
    output:
        raw = "results/magma/output/{pheno}/{geno_set}/chr{chr}.genes.raw",
        out = "results/magma/output/{pheno}/{geno_set}/chr{chr}.genes.out"
    params:
        prefix_out = lambda x, output: Path(str(output[0])).with_suffix('').with_suffix('').as_posix(),
        prefix_bfile = lambda x, input: Path(str(input.bfile[0])).with_suffix('').as_posix(),
        n = lambda x: get_n_trait(x)
    shell:
        """
        magma \
	        --bfile {params.prefix_bfile} \
	        --gene-annot {input.genes_annot} \
	        --pval {input.gwas_file} use=SNP,p N={params.n} \
	        --gene-model snp-wise=mean \
	        --out {params.prefix_out}
        """
    
# aggregate results
rule magma_agg_gene_assoc:
    input:
        out = expand("results/magma/output/{{pheno}}/{{geno_set}}/chr{chr}.genes.out", \
               chr = range(1,23)),
        raw = expand("results/magma/output/{{pheno}}/{{geno_set}}/chr{chr}.genes.raw", \
               chr = range(1,23)),
    output:
        raw = "results/magma/aggregate/{geno_set}/{pheno}.genes.raw",
        out = "results/magma/aggregate/{geno_set}/{pheno}.genes.out"
    shell:
        """
        awk 'NR == FNR {{print; next}} FNR > 1' {input.out} > {output.out}
        awk 'NR == FNR {{print; next}} FNR > 2' {input.raw} > {output.raw}
		"""

# update gene out with hgnc_symbols & fdr correction (with qvalue)
rule magma_format_output:
    input:
        out = rules.magma_agg_gene_assoc.output.out,
        gene_loc = "resources/magma/geneanno.loc"
    output:
        "results/magma/aggregate/{geno_set}/{pheno}.magma.tsv"
    conda: "renv"
    script:
        "snakescripts/magma_format_output.R"

# run results for all trait
def get_trait_input(x):
    """
    Get the input for all traits
    """
    df_meta = pd.read_csv("resources/metadata/trait_metadata_n.tsv", encoding_errors = "ignore", sep="\t")
    pheno_meta = df_meta['trait_id'].unique()
    input_dir = Path("resources/ma/")
    pheno = [f.stem for f in input_dir.glob("*.ma") \
             if f.is_file() and f.stem in pheno_meta]
    files = [f"results/magma/aggregate/{x.geno_set}/{p}.magma.tsv" for p in pheno]
    return files

rule magma_all_trait:
    input: get_trait_input
    output: "results/aggregate/{geno_set}.magma.gz.parquet"
    conda: "renv"
    log: "logs/aggregate/{geno_set}.magma.log"
    script: "snakescripts/aggregate/magma.R"