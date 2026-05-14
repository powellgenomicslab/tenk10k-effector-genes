# run Genetic correlations analysis to get N independent traits

rule ldak_calc_tag_file:
    input:
        bfiles = expand("resources/genotypes/{{study}}/chr{{chr}}.{ext}",
                        ext = ["bed", "bim", "fam"], chr = range(1, 23)),
    output:
        "resources/ldak/tagfiles/{study}/chr{chr}.tagging"
    params:
        bfile_prefix = "resources/genotypes/{study}/chr{chr}",
        power = "-.25",
        out_prefix = "resources/ldak/tagfiles/{study}/chr{chr}"
    log: "logs/ldak/{study}/chr{chr}.tag.log"
    shell:
        """
        ldak \
            --calc-tagging {params.out_prefix} \
            --bfile {params.bfile_prefix} \
            --power {params.power}
        """

rule ldak_join_tag_files:
    input:
        tags = expand("resources/ldak/tagfiles/{{study}}/chr{chr}.tagging",
                      chr = range(1, 23)),
    output:
        taglist = "resources/ldak/tag_merged/{study}.taglist.txt",
        tagging = "resources/ldak/tag_merged/{study}.tagging"
    log: "logs/ldak/{study}/merge_tag.log"
    params:
        out_prefix = "resources/ldak/tag_merged/{study}"
    shell:
        """
        mkdir -p resources/ldak/tag_merged
        echo {input.tags} | tr ' ' '\\n' > {output.taglist}
        ldak --join-tagging {params.out_prefix} --taglist {output.taglist}
        """

rule ldak_format_ma_sumstats:
    input:
        ma = "resources/ma/{phenotype}.ma",
        trait_metadata = "resources/metadata/trait_metadata_curated.xlsx"
    output:
        "resources/ldak/sumstats/{phenotype}.sumstats.ldak"
    log: "logs/ldak/format_sumstats/{phenotype}.ma"
    conda: "renv"
    script: "snakescripts/ldak/format_ma_sumstats.R"

def all_traits_sumstats_ldak():
    files = list(Path("resources/ma/").glob("*.ma"))
    targets = [f"resources/ldak/sumstats/{f.stem}.sumstats.ldak" \
               for f in files if f.is_file()]
    return targets

rule all_traits_sumstats_ldak:
    input: all_traits_sumstats_ldak()
    output: touch("resources/ldak/.done/all_sumstats.done")
    shell: "touch {output}"

rule ldak_gen_cor:
    """
    Calculate genetic correlation using LDAK
    """
    input:
        sumstats1 = "resources/ldak/sumstats/{trait1}.sumstats.ldak",
        sumstats2 = "resources/ldak/sumstats/{trait2}.sumstats.ldak",
        tag_file = "resources/ldak/tag_merged/{study}.tagging",
    output:
        "results/gen_cor/{study}/{trait1}.{trait2}.cors"
    params:
        out_prefix = "results/gen_cor/{study}/{trait1}.{trait2}",
        allow_ambiguous = "YES",
        check_sums = "NO"
    threads: 8
    shell:
        """
        mkdir -p results/gen_cor/{wildcards.study}
        ldak --sum-cors {params.out_prefix} \
            --summary {input.sumstats1} \
            --summary2 {input.sumstats2} \
            --tagfile {input.tag_file} \
            --max-threads {threads} \
            --check-sums {params.check_sums} \
            --allow-ambiguous {params.allow_ambiguous}
        """


rule mk_trait_list:
    output: "resources/ldak/trait_list/{study}.txt"
    conda: "renv"
    script: "snakescripts/prep_trait_list/{wildcards.study}.R"

rule mk_pairwise_trait:
    """
    make pairwise combination of traits
    """
    input: "resources/ldak/trait_list/{study}.txt"
    output: "resources/ldak/misc/{study}.trait_pairs.gen_cor.txt"
    params:
        ldak_prefix = "results/gen_cor/{study}"
    conda: "renv"
    script: "snakescripts/ldak/make_pairwise_trait_gencor.R"

def ldak_gen_cor_all(x):
    with open(f"resources/ldak/misc/{x.study}.trait_pairs.gen_cor.txt") as f:
        files = [line.strip() for line in f]
    return files


rule all_gencor_ldak:
    input: ldak_gen_cor_all
    output: "results/aggregate/{study}.gen_cor.ldak.tsv"
    conda: "renv"
    log: "logs/ldak/gen_cor/{study}.all.log"
    script: "snakescripts/aggregate/gen_cor_ldak.R"

rule gencor_eigen:
    """create eigenvectors & eigen values for genetic correlation"""
    input: "results/aggregate/{study}.gen_cor.ldak.tsv"
    output: "results/aggregate/{study}.gen_cor_eigen.tsv"
    conda: "renv"
    log: "logs/gen_cor/{study}.eigen.log"
    script: "snakescripts/aggregate/gen_cor_eigen.R"

# rule ldak_gen_cor_all:
#     input: ldak_gen_cor_all
#     output: touch("results/.done/gen_cor.{study}.done")
#     shell: "mkdir -p $(dirname {output}) && touch {output}"

# rule ldsc_munge_sumstats:
#     input:
#         script = "softwares/ldsc/munge_sumstats.py"
#     output:
#         "results/gen_cor/{trait}_munge_sumstats.txt"
#     params:
#         trait="{trait}"
#     log:
#         "logs/gen_cor/{trait}.log"
#     shell:
#         """
#         echo "Running munge sumstats for {params.trait}" > {log}
#         # Simulate the command for munge sumstats
#         echo "Munged sumstats for {params.trait}" > {output}
#         echo "Munge completed for {params.trait}" >> {log}
#         """

# rule ldsc_rho:
#     input:
#         script = "softwares/ldsc/ldsc.py"
#     output:
#         "results/gen_cor/{trait}_gen_cor.txt"
#     params:
#         trait="{trait}"
#     log:
#         "logs/gen_cor/{trait}.log"
#     shell:
#         """
#         echo "Running genetic correlation analysis for {params.trait}" > {log}
#         # Simulate the command for genetic correlation analysis
#         echo "Genetic correlation results for {params.trait}" > {output}
#         echo "Analysis completed for {params.trait}" >> {log}
#         """