rule run_ivw:
    input:
        pheno_file = config["workflow"]["phenotype_dir"] + "/{trait}.txt"
    output:
        temp(touch(config["workflow"]["ivw_output_dir"] + "/.{trait}.done"))
    params:
        ivw_output = config["workflow"]["ivw_output_dir"],
        smr = config["container_paths"]["sensitivity_smr_results"],
        snp_smr = config["container_paths"]["snp_smr_results"]
    resources:
        ncpus = 8,
        mem = "64G",
        time = "23:00:00",
        jobfs = "1G",
        queue = "express"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        while IFS= read -r CELL_TYPE || [ -n "$CELL_TYPE" ]; do
            [ -z "$CELL_TYPE" ] && continue
            mkdir -p {params.ivw_output}/output/$CELL_TYPE/
            Rscript sensitivity/rules/snakescripts/ivw_mr/run_ivw_corr.R \
                {wildcards.trait} "$CELL_TYPE" \
                {params.ivw_output}/output/$CELL_TYPE/{wildcards.trait}_IVW-LD_results.tsv \
                {params.smr} {params.snp_smr} \
                || {{ echo "Skipping {wildcards.trait} - $CELL_TYPE"; continue; }}
        done < {input.pheno_file}
        """


rule process_ivw:
    input:
        config["workflow"]["ivw_output_dir"] + "/.{trait}.done"
    output:
        temp(touch(config["workflow"]["ivw_output_dir"] + "/.{trait}.processed"))
    params:
        ivw_output = config["workflow"]["ivw_output_dir"]
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        "mkdir -p {params.ivw_output}/results && Rscript sensitivity/rules/snakescripts/ivw_mr/process_ivw_corr.R {wildcards.trait}"


def get_all_processed(wildcards):
    checkpoints.prepare_phenotype_list.get()
    with open(config["workflow"]["phenotype_list"]) as f:
        traits = [line.strip() for line in f if line.strip()]
    return expand(
        config["workflow"]["ivw_output_dir"] + "/.{trait}.processed",
        trait=traits
    )


rule make_ivw_parquet:
    input:
        get_all_processed
    output:
        touch(config["workflow"]["ivw_output_dir"] + "/.make_parquet.done")
    params:
        ivw_output = config["workflow"]["ivw_output_dir"]
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        Rscript sensitivity/rules/snakescripts/ivw_mr/make_ivw_corr_parquet.R
        rm -rf {params.ivw_output}/output/
        rm -f {params.ivw_output}/results/*_IVW-LD_results.tsv
        """
