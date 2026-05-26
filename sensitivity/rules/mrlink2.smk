def get_mrlink2_exposure_sentinels(wildcards):
    checkpoints.prepare_phenotype_list.get()
    pheno_file = config["workflow"]["phenotype_dir"] + f"/{wildcards.trait}.txt"
    with open(pheno_file) as f:
        cell_types = [line.strip() for line in f if line.strip()]
    return expand(
        config["workflow"]["mrlink2_output_dir"] + "/inputs/exposures/eqtl/.{cell_type}.done",
        cell_type=cell_types
    )


def get_all_mrlink2_processed(wildcards):
    checkpoints.prepare_phenotype_list.get()
    with open(config["workflow"]["phenotype_list"]) as f:
        traits = [line.strip() for line in f if line.strip()]
    return expand(
        config["workflow"]["mrlink2_output_dir"] + "/outputs/.{trait}.processed",
        trait=traits
    )


rule preprocess_mrlink2_outcome:
    input:
        config["workflow"]["phenotype_dir"] + "/{trait}.txt"
    output:
        temp(touch(config["workflow"]["mrlink2_output_dir"] + "/inputs/outcomes/.{trait}.done"))
    params:
        mrlink2_output = config["workflow"]["mrlink2_output_dir"],
        ma_dir = config["container_paths"]["ma_dir"]
    resources:
        ncpus = 4,
        mem = "32G",
        time = "23:00:00",
        jobfs = "1G",
        queue = "express"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        mkdir -p {params.mrlink2_output}/inputs/outcomes
        Rscript sensitivity/rules/snakescripts/mrlink2/preprocess_mrlink2_outcome.R \
            {wildcards.trait} \
            {params.ma_dir}/{wildcards.trait}.ma \
            {params.mrlink2_output}/inputs/outcomes/{wildcards.trait}.parq
        """


rule preprocess_mrlink2_exposure:
    input:
        config["key_files"]["smr_results"]
    output:
        temp(touch(config["workflow"]["mrlink2_output_dir"] + "/inputs/exposures/eqtl/.{cell_type}.done"))
    params:
        mrlink2_output = config["workflow"]["mrlink2_output_dir"],
        smr = config["container_paths"]["smr_results"]
    resources:
        ncpus = 4,
        mem = "32G",
        time = "06:00:00",
        jobfs = "1G",
        queue = "normal"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        mkdir -p {params.mrlink2_output}/inputs/exposures/eqtl/{wildcards.cell_type}
        Rscript sensitivity/rules/snakescripts/mrlink2/preprocess_mrlink2_exposure.R \
            {params.smr} \
            {wildcards.cell_type} \
            {params.mrlink2_output}/inputs/exposures/eqtl/{wildcards.cell_type}
        """


rule preprocess_mrlink2_exposure_traits:
    input:
        pheno_file = config["workflow"]["phenotype_dir"] + "/{trait}.txt"
    output:
        temp(touch(config["workflow"]["mrlink2_output_dir"] + "/inputs/exposures/phenotypes/.{trait}.done"))
    params:
        mrlink2_output = config["workflow"]["mrlink2_output_dir"],
        smr = config["container_paths"]["smr_results"]
    resources:
        ncpus = 4,
        mem = "32G",
        time = "02:00:00",
        jobfs = "1G",
        queue = "normal"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        mkdir -p {params.mrlink2_output}/inputs/exposures/phenotypes/{wildcards.trait}
        while IFS= read -r CELL_TYPE || [ -n "$CELL_TYPE" ]; do
            [ -z "$CELL_TYPE" ] && continue
            Rscript sensitivity/rules/snakescripts/mrlink2/preprocess_mrlink2_exposure_traits.R \
                {wildcards.trait} "$CELL_TYPE" \
                {params.smr} \
                {params.mrlink2_output}/inputs/exposures/phenotypes/{wildcards.trait}/${{CELL_TYPE}}_mrlink2_genes.tsv
        done < {input.pheno_file}
        """


rule run_mrlink2:
    input:
        pheno_file = config["workflow"]["phenotype_dir"] + "/{trait}.txt",
        outcome = config["workflow"]["mrlink2_output_dir"] + "/inputs/outcomes/.{trait}.done",
        exposure_traits = config["workflow"]["mrlink2_output_dir"] + "/inputs/exposures/phenotypes/.{trait}.done",
        exposures = get_mrlink2_exposure_sentinels
    output:
        temp(touch(config["workflow"]["mrlink2_output_dir"] + "/outputs/.{trait}.done"))
    params:
        mrlink2_output = config["workflow"]["mrlink2_output_dir"]
    resources:
        ncpus = 8,
        mem = "32G",
        time = "23:00:00",
        jobfs = "1G",
        queue = "express"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        while IFS= read -r CELL_TYPE || [ -n "$CELL_TYPE" ]; do
            [ -z "$CELL_TYPE" ] && continue
            GENES_FILE={params.mrlink2_output}/inputs/exposures/phenotypes/{wildcards.trait}/${{CELL_TYPE}}_mrlink2_genes.tsv
            [ ! -f "$GENES_FILE" ] && {{ echo "Skipping {wildcards.trait} - $CELL_TYPE: no genes file"; continue; }}
            while IFS=$'\t' read -r CHROM probeID || [ -n "$CHROM" ]; do
                [ -z "$CHROM" ] && continue
                OUTPUT_DIR={params.mrlink2_output}/outputs/$CELL_TYPE/{wildcards.trait}/chr$CHROM
                mkdir -p $OUTPUT_DIR
                if [ -f "$OUTPUT_DIR/$probeID.txt" ] || [ -f "$OUTPUT_DIR/$probeID.txt_noregions" ]; then
                    echo "Skipping $probeID chr$CHROM - already exists"
                    continue
                fi
                python3 /workspace/mrlink2/mr_link_2_standalone.py \
                    --reference_bed /genotypes/chr${{CHROM}}_common_variants \
                    --sumstats_exposure {params.mrlink2_output}/inputs/exposures/eqtl/$CELL_TYPE/chr$CHROM/$probeID.txt \
                    --sumstats_outcome {params.mrlink2_output}/inputs/outcomes/{wildcards.trait}.parq \
                    --out $OUTPUT_DIR/$probeID.txt \
                    --maf_threshold 0.01 \
                    --no_exclude_hla \
                    || echo "Skipping $probeID chr$CHROM - $CELL_TYPE: mr_link_2 failed"
            done < <(tail -n +2 "$GENES_FILE")
        done < {input.pheno_file}
        """


rule process_mrlink2:
    input:
        pheno_file = config["workflow"]["phenotype_dir"] + "/{trait}.txt",
        run_done = config["workflow"]["mrlink2_output_dir"] + "/outputs/.{trait}.done"
    output:
        temp(touch(config["workflow"]["mrlink2_output_dir"] + "/outputs/.{trait}.processed"))
    params:
        mrlink2_output = config["workflow"]["mrlink2_output_dir"]
    resources:
        ncpus = 4,
        mem = "32G",
        time = "00:10:00",
        jobfs = "1G",
        queue = "express"
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        while IFS= read -r CELL_TYPE || [ -n "$CELL_TYPE" ]; do
            [ -z "$CELL_TYPE" ] && continue
            INPUT_DIR={params.mrlink2_output}/outputs/$CELL_TYPE/{wildcards.trait}
            OUTPUT_FILEPATH={params.mrlink2_output}/outputs/$CELL_TYPE/{wildcards.trait}_mrlink2_aggregated.tsv
            Rscript sensitivity/rules/snakescripts/mrlink2/process_mrlink2_results.R \
                {wildcards.trait} "$CELL_TYPE" \
                "$INPUT_DIR" "$OUTPUT_FILEPATH" \
                || echo "Skipping {wildcards.trait} - $CELL_TYPE: process failed"
        done < {input.pheno_file}
        """


rule combine_mrlink2_results:
    input:
        get_all_mrlink2_processed
    output:
        touch(config["workflow"]["mrlink2_output_dir"] + "/.combine.done")
    params:
        mrlink2_output = config["workflow"]["mrlink2_output_dir"]
    singularity:
        config.get("mr_singularity_image", "")
    shell:
        """
        mkdir -p {params.mrlink2_output}/results
        Rscript sensitivity/rules/snakescripts/mrlink2/combine_mrlink2_results.R \
            {params.mrlink2_output}/outputs \
            {params.mrlink2_output}/results/tenk10k_phase1_sensitivity.mrlink2.parquet.gz
        find {params.mrlink2_output}/outputs -mindepth 2 -maxdepth 2 -type d -exec rm -rf {{}} \;
        rm -rf {params.mrlink2_output}/inputs/exposures/
        """
