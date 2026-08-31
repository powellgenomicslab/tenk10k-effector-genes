# snakemake rules to format GWAS summary stats

# NOTE ON GWAS PROVENANCE
# `resources/ma/` is a curated *input* to this pipeline, not one of its outputs.
# Every downstream arm (SMR, MAGMA, coloc, LDAK, and the sensitivity MR pipeline)
# reads `resources/ma/{pheno}.ma`, and no rule writes it. The directory is
# assembled by hand from three routes:
#   1. `rule format_gwas` below, for the 25 traits with a trait-specific script
#      in snakescripts/format_gwas/. It writes to `resources/pipeline_ma/`, a
#      staging area, and the harmonised files are reviewed before being promoted
#      into `resources/ma/`.
#   2. `rule extract_finngen_gwas` + `rule format_finngen_gwas`, for the
#      phenotypes listed in workflow/config/finngen_meta_path.yaml.
#   3. Summary statistics already supplied in GCTA .ma format by the source
#      consortium, copied in directly.
# Anyone reproducing this pipeline must assemble `resources/ma/` the same way;
# see the "GWAS inputs" section of workflow/README.md.

rule format_gwas:
    """
    Format GWAS summary statistics into GCTA .ma format.

    Writes to the resources/pipeline_ma/ staging area; see the provenance note
    above for how those files reach resources/ma/.
    """
    input:
        gwas = "resources/sumstats/gwas/{pheno}.gwas",
        trait_metadata = "metadata/trait.tsv",
        hg19tohg38 = "resources/misc/hg19ToHg38.over.chain",
        liftover_script = "workflow/rules/snakescripts/hg19tohg38.R"
    output: "resources/pipeline_ma/{pheno}.ma"
    threads: 8
    resources:
        mem = "16G",
        jobfs = "8G",
        ncpus = 8
    conda: "../envs/renv.yaml"
    log: "logs/format_gwas/{pheno}.log"
    script: "snakescripts/format_gwas/{wildcards.pheno}.R"

rule mk_tabix_region:
    """
    Make tabix region file based on genotype file for extraction
    """
    input: expand("resources/genotypes/{{study}}_common/chr{chr}.bim", chr = range(1,23))
    output: "resources/misc/tabix_region/{study}_common.region"
    threads: 8
    shell: "cut -f1,4 {input} | sort -nk1 -nk2 > {output}"

def get_finngen_gwas(x):
    with open("workflow/config/finngen_meta_path.yaml", 'r') as f:
        config = yaml.full_load(f)
    base_dir = Path(config['base_dir'])
    return {
        "gwas": str(base_dir / f"{config['pheno_prefix'][x.pheno]}_meta_out.tsv.gz"),
        "tbi": str(base_dir / f"{config['pheno_prefix'][x.pheno]}_meta_out.tsv.gz.tbi")
    }

rule extract_finngen_gwas:
    input:
        unpack(get_finngen_gwas),
        region = "resources/misc/tabix_region/{study}_common.region"
    output: temp("resources/sumstats/finngen_gwas_extract/{study}/{pheno}.tsv.gz")
    threads: 8
    envmodules: "htslib"
    params:
        # extract chr pos ref alt snp, af cohorts, estimate meta (all), rsid
        # cols = ",".join([str(x) for x in [1, 2, 3, 4, 5, 9, 16, 20, 25, 30, 31, 32, 33, 34, 67]])
    shell:
        """
        tabix -hR {input.region} {input.gwas} | \
            gzip -c > {output}
        """

# rule mk_tabix_region:
#     """
#     Make tabix region file based on genotype file for extraction
#     """
#     input: expand("resources/genotypes/{{study}}_common/chr{chr}.bim", chr = range(1,23))
#     output: "resources/misc/tabix_region/{{study}}_common.region"
#     threads: 8
#     envmodules: "htslib"
#     shell:

rule format_finngen_gwas:
    """
    Format GWAS summary statistics from Finggen meta-analysis (with UKB + MVP EUR)
    for use in SMR
    """
    input:
        unpack(get_finngen_gwas),
        gwas_subset = rules.extract_finngen_gwas.output,
        metadata = "resources/misc/meta_analysis_mvp_ukbb_FinnGen_R12_MVP_UKBB_mapping.tsv",
        frq = "resources/genotypes_frq/{study}.frq",
    output: "resources/sumstats/finngen_gwas_extract/{study}/{pheno}.ma"
    threads: 8
    conda: "../envs/renv.yaml"
    log: "logs/format_finngen_gwas/{study}/{pheno}.log"
    script: "snakescripts/format_finngen_gwas.R"

rule split_ma_by_chr:
    """
    Split GWAS ma format by chromosome
    """
    input: "resources/ma/{pheno}.ma"
    output: expand("resources/ma_by_chr/{{pheno}}/chr{chr}.ma", chr = range(1,23))
    shell:
        """
        OUTDIR=$(dirname {output[0]})
        mkdir -p $OUTDIR && \
            awk -F'\\t' -v OFS='\\t' -v OUTDIR="$OUTDIR" '
                NR==1 {{ $7='P'; header=$0; next }}
                {{split($1, a, ":"); 
                 chr=a[1]; 
                 file=OUTDIR "/chr" chr ".ma"; 
                 if (!seen[chr]++) {{print header > file;}}
                 print $0 >> file;
                }}' {input}
        """