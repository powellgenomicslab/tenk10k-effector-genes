import os

if 'PBS_JOBID' in os.environ or 'SLURM_JOB_ID' in os.environ:
    configfile: "sensitivity/config/gadi.yaml"
else:
    configfile: "sensitivity/config/local.yaml"

include: "rules/ivw_mr.smk"
include: "rules/mrlink2.smk"


checkpoint prepare_phenotype_list:
    input:
        config["key_files"]["smr_results"]
    output:
        pheno_dir = directory(config["workflow"]["phenotype_dir"]),
        sentinel = config["workflow"]["phenotype_list"]
    params:
        smr = config["container_paths"]["smr_results"],
        output_path = config["workflow"]["phenotype_dir"] + "/"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        mkdir -p {output.pheno_dir}
        Rscript sensitivity/scripts/preparePhenotypeList.R \
            {params.smr} {params.output_path} {output.sentinel}
        """

rule all_ivw:
    input:
        config["workflow"]["ivw_output_dir"] + "/.make_parquet.done"


rule all_mrlink2:
    input:
        config["workflow"]["mrlink2_output_dir"] + "/.combine.done"
