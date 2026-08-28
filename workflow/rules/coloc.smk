# Rules for coloc & multivariant coloc analysis

import math
from pathlib import Path

configfile: "workflow/config/mvcoloc.yaml"
# ======================================================================
# Standard coloc (per biosample × phenotype × chromosome)
# ======================================================================

checkpoint prep_coloc_input:
    """Prepare coloc input files (eQTL data)"""
    output: directory("resources/coloc/{study}/")
    conda: "../envs/renv.yaml"
    script: "snakescripts/prep_coloc_input/{wildcards.study}.R"


rule run_coloc:
    """Run coloc per chr / cell type / phenotype"""
    input:
        eqtl = ancient("resources/coloc/{study}/{biosample}/common_egenes/chr{chr}.fst"),
        gwas = "resources/ma_by_chr/{pheno}/chr{chr}.ma",
        gene_loc = ancient("resources/misc/gencode.v44.gene_type.tsv"),
        pheno_metadata = ancient("metadata/trait.tsv")
    output:
        coloc = temp("results/coloc/{study}/{biosample}/{pheno}/chr{chr}.coloc.tsv")
    # Serial workload: the per-gene loop is map_df(egene_list, run_coloc) and
    # coloc.abf is single-threaded, so extra cores sit idle. Measured on the
    # ASDC/sle pilot: 16 CPU-seconds against 39 wall-seconds on 4 requested
    # cores (<0.5 core used). Requesting 4 cores charged ~4x the SU for nothing.
    # Peak memory was 2.9 GB, so 8G is kept as headroom for larger biosamples.
    threads: 1
    log: "logs/coloc/{study}/{biosample}/{pheno}/chr{chr}.log"
    resources:
        mem = "8G",
        ncpus = 1
    params:
        window_bp = 100000
    conda: "../envs/renv.yaml"
    script: "snakescripts/coloc/run_coloc_test.R"


rule concat_coloc_chr:
    """Concatenate coloc results across chromosomes"""
    input:
        coloc = expand(
            "results/coloc/{{study}}/{{biosample}}/{{pheno}}/chr{chr}.coloc.tsv",
            chr=range(1, 23)
        )
    output:
        coloc = "results/coloc/{study}/{biosample}/{pheno}/all_chr.coloc.tsv"
    shell:
        # Header-aware concat. The previous `awk 'NR == 1 || FNR > 1'` emitted
        # only the first chromosome's header and then appended rows of any
        # width, so files written by different versions of run_coloc_test.R were
        # silently interleaved (12- vs 13-column layouts) and the SNP count
        # ended up in an unnamed field. Fail loudly instead.
        # NB: literal awk braces must be doubled - snakemake runs the shell
        # string through str.format() (same as split_ma_by_chr in format_gwas.smk).
        r"""
        awk -F'\t' '
            FNR == 1 {{
                if (NR == 1) {{ header = $0; print }}
                else if ($0 != header) {{
                    printf "ERROR: header mismatch in %s\n  expected: %s\n  found:    %s\n", \
                        FILENAME, header, $0 > "/dev/stderr"
                    exit 1
                }}
                next
            }}
            {{ print }}
            END {{
                if (header == "") {{
                    print "ERROR: no header found - all inputs were empty" > "/dev/stderr"
                    exit 1
                }}
            }}
        ' {input.coloc} > {output.coloc}
        """


def target_coloc(x, coloc_ext):
    """Target function for all biosamples × phenotypes"""
    dir_study = checkpoints.prep_coloc_input.get(**x).output[0]
    BIOSAMPLES = [d.name for d in Path(dir_study).iterdir() if d.is_dir()]
    with open("resources/misc/target_phenotypes.txt") as f:
        PHENOS = [line.strip() for line in f if line.strip()]
    return expand(
        f"results/coloc/{x.study}/{{biosample}}/{{pheno}}/all_chr.{coloc_ext}.tsv",
        biosample=BIOSAMPLES, pheno=PHENOS
    )


rule concat_coloc_all:
    """Aggregate coloc results for all biosamples and phenotypes"""
    input: lambda x: target_coloc(x, "coloc")
    # NOTE: this rule produces the `v2` aggregate. The published summary
    # statistics are built from `v3`, which is `v2` plus coloc recomputed for the
    # 34 traits whose GWAS carried zero-SE or duplicate-ID SNPs (11,930,730 rows
    # -> 15,434,261 rows; see CHANGELOG.v5.md). That recomputation was run out of
    # band and has no rule here, so re-running this arm reproduces v2, not v3.
    output: "results/aggregate/coloc/{study}.coloc.v2.parquet.gz"
    params:
        schema_helper = "workflow/rules/snakescripts/aggregate/coloc_schema.R"
    conda: "../envs/renv.yaml"
    log: "logs/aggregate/{study}.concat_coloc.log"
    resources:
        mem = "64G",
        ncpus = 8
    script: "snakescripts/aggregate/concat_coloc_all.R"


# ======================================================================
# LD matrices for multivariant coloc
# ======================================================================

rule make_ld_eqtl:
    """Make LD matrix per chromosome for eQTL data"""
    input:
        bfile = expand(
            "resources/genotypes/{{study}}/chr{{chr}}.{ext}",
            ext=["bed", "bim", "fam"]
        ),
        eqtl_dir = ancient("resources/coloc/{study}"),
        gene_loc = ancient("resources/misc/gencode.v44.gene_type.tsv")
    output:
        ld = directory("resources/ld/eqtl/{study}/chr{chr}")
    threads: 8
    log: "logs/ld/{study}/chr{chr}.log"
    resources:
        mem = "32G",
        ncpus = 8
    params:
        window_bp = 100000
    conda: "../envs/renv.yaml"
    script: "snakescripts/make_ld/eqtl/{wildcards.study}.R"


# Run SuSiE GWAS by chromosome for multivariant coloc
def _get_traits():
    with open("resources/misc/target_phenotypes.txt") as f:
        traits = [line.strip() for line in f if line.strip()]
    return traits

TRAITS = _get_traits()

checkpoint generate_susie_gwas_cmds:
    """
    Generate per-task shell scripts and a command file listing them.
    """
    input:
        gene_loc = ancient("resources/misc/gencode.v44.gene_type.tsv"),
        pheno_metadata = ancient("metadata/trait.tsv"),
        dir_bfile = ancient("resources/genotypes/{study}"),
        script = "workflow/rules/snakescripts/coloc/run_susie_gwas_cli.R",
        gwas = expand("resources/ma_by_chr/{t}/chr{c}.ma", t=TRAITS, c=range(1, 23))
    output:
        cmd_file = "cmds/susie_gwas/{study}/main_cmds.txt"
    params:
        traits = lambda wc: " ".join(TRAITS),
        outdir = "results/susie_gwas/{study}",
        window_bp = 100000,
        runsusie_coverage = 0.1,
        min_p_gwas = 1e-4,
        runsusie_maxit = 200,
        runsusie_repeat = False,
        nthreads = 1
    conda: "../envs/renv.yaml"
    shell:
        r"""
        PROJECT_ROOT=$(pwd)
        CONDA_ENV="$CONDA_PREFIX"
        SCRIPT_DIR=$(dirname {output.cmd_file})

        mkdir -p "$SCRIPT_DIR" > {output.cmd_file}

        for chr in $(seq 1 22); do
            for pheno in {params.traits}; do
                outfile="{params.outdir}/${{pheno}}/chr${{chr}}.susie.rds"
                script_file="$SCRIPT_DIR/task_chr${{chr}}_${{pheno}}.sh"

                cat > "$script_file" <<TASK_EOF
#!/bin/bash
set -euo pipefail
cd $PROJECT_ROOT
if [ -f "$outfile" ] && [ -s "$outfile" ]; then
    echo "\$(date): SKIP $outfile"
    exit 0
fi
export PATH=$CONDA_ENV/bin:\$PATH
export CONDA_PREFIX=$CONDA_ENV
export R_LIBS_SITE=$CONDA_ENV/lib/R/library

mkdir -p \$(dirname $outfile)

Rscript {input.script} \
    --study {wildcards.study} \
    --chr $chr \
    --pheno $pheno \
    --gwas resources/ma_by_chr/${{pheno}}/chr${{chr}}.ma \
    --gene_loc {input.gene_loc} \
    --pheno_metadata {input.pheno_metadata} \
    --dir_bfile {input.dir_bfile} \
    --output $outfile \
    --threads {params.nthreads} \
    --window_bp {params.window_bp} \
    --runsusie_coverage {params.runsusie_coverage} \
    --min_p_gwas {params.min_p_gwas} \
    --runsusie_maxit {params.runsusie_maxit}
TASK_EOF

                chmod +x "$script_file"
                echo "$PROJECT_ROOT/$script_file" >> {output.cmd_file}
            done
        done

        echo "Generated $(wc -l < {output.cmd_file}) task scripts"
        """

rule run_susie_gwas:
    """
    Execute SuSiE GWAS tasks via nci-parallel.
    """
    input:
        cmd_file = "cmds/susie_gwas/{study}/main_cmds.txt"
    output:
        done = touch("results/susie_gwas/{study}/.done")
    resources:
        ncpus  = 96,
        mem    = "380G",
        time   = "24:00:00",
        jobfs  = "20G"
    params:
        ncores_per_task        = 1,  
        ncores_per_numanode     = 12, # 48  → ppr
        timeout         = 7200
    log: "logs/susie_gwas/{study}/nciparallel.log"
    shell:
        """
        module load nci-parallel/1.0.0a 2>> {log}

        NCMDS=$(wc -l < {input.cmd_file})

        M=$(( PBS_NCI_NCPUS_PER_NODE / {params.ncores_per_task} )) 
        PPR=$(( {params.ncores_per_numanode} / {params.ncores_per_task} ))
        NP=$((M * PBS_NCPUS / PBS_NCI_NCPUS_PER_NODE))
        
        echo $PPR
        echo $NP
        mpirun -np $NP \
            --map-by ppr:${{PPR}}:NUMA:PE={params.ncores_per_task} \
            nci-parallel \
            --verbose \
            --input-file {input.cmd_file} \
            --timeout {params.timeout} \
        >> {log} 2>&1 \
        && touch {output.done}
        """

rule rerun_susie_gwas:
    input:
        dir_res = "results/susie_gwas/{study}",
        main_cmds = "cmds/susie_gwas/{study}/main_cmds.txt"
    output:
        rerun_cmds = "cmds/susie_gwas/{study}/rerun_cmds.txt",
        done = touch("results/susie_gwas/{study}/.done.rerun")
    log: "logs/susie_gwas/{study}/nciparallel.rerun.log"
    resources:
        ncpus  = 96,
        mem    = "380G",
        time   = "24:00:00",
        jobfs  = "20G"
    params:
        traits = lambda wc: " ".join(TRAITS),
        outdir = "results/susie_gwas/{study}",
        ncores_per_task        = 2,  
        ncores_per_numanode     = 12, # 48  → ppr
        timeout         = 7200
    conda: "../envs/renv.yaml"
    shell:
        """
        module load nci-parallel/1.0.0a 2>> {log}
        PROJECT_ROOT=$(pwd)
        CONDA_ENV="$CONDA_PREFIX"
        SCRIPT_DIR=$(dirname {output.rerun_cmds})
        mkdir -p "$SCRIPT_DIR" > {output.rerun_cmds}
        
        for chr in $(seq 1 22); do
            for pheno in {params.traits}; do
                outfile="{params.outdir}/${{pheno}}/chr${{chr}}.susie.rds"
                if [ -f "$outfile" ] && [ -s "$outfile" ]; then
                    continue
                fi
                script_file="$SCRIPT_DIR/task_chr${{chr}}_${{pheno}}.sh"
                echo "$PROJECT_ROOT/$script_file" >> {output.rerun_cmds}
            done
        done

        echo "Generated $(wc -l < {output.rerun_cmds}) task scripts to rerun" 2>&1 | tee -a {log}

        M=$(( PBS_NCI_NCPUS_PER_NODE / {params.ncores_per_task} ))
        PPR=$(( {params.ncores_per_numanode} / {params.ncores_per_task} ))
        NP=$((M * PBS_NCPUS / PBS_NCI_NCPUS_PER_NODE))
        
        mpirun -np $NP \
            --map-by ppr:${{PPR}}:NUMA:PE={params.ncores_per_task} \
            nci-parallel \
            --verbose \
            --input-file {output.rerun_cmds} \
            --timeout {params.timeout} \
        2>&1 | tee -a {log} \
        && touch {output.done}
        """

rule susie_gwas_rds_check:
    """
    Bridge: individual RDS depends on the batch sentinel from run_susie_gwas.
    Checks that the RDS was actually created and is a valid R object.
    Creates a header-only placeholder if no GWAS signal for this chr×trait.
    """
    input:
        done = "results/susie_gwas/{study}/.done"  
    output:
        rds = "results/susie_gwas/{study}/{pheno}/chr{chr}.susie.rds"
    conda: "../envs/renv.yaml"
    shell:
        """
        # Force NFS cache refresh before checking
        ls -la $(dirname {output.rds}) > /dev/null 2>&1

        if [ ! -f {output.rds} ]; then
            echo "ERROR: RDS missing after batch completed: {output.rds}" >&2
            echo "  Check logs for chr={wildcards.chr} pheno={wildcards.pheno}" >&2
            exit 1
        fi

        # Validate it is a readable RDS (not a truncated/corrupt file)
        Rscript -e "
            rds <- tryCatch(readRDS('{output.rds}'), error = function(e) stop(e))
            if (!is.list(rds)) stop('RDS is not a list: ', class(rds))
            n_genes <- sum(!sapply(rds, is.null))
            cat(sprintf('OK: {wildcards.pheno} chr{wildcards.chr} — %d/%d genes with SuSiE results\n',
                        n_genes, length(rds)))
        " || {{ echo "ERROR: RDS is corrupt or unreadable: {output.rds}" >&2; exit 1; }}

        touch {output.rds}
        """

rule generate_mvcoloc_cmds:
    """
    Generate command scripts for multivariate coloc analysis.
    """
    input:
        susie_done = "results/susie_gwas/{study}/.done",
        dir_ld_eqtl = "resources/ld/eqtl/{study}",
        dir_eqtl = "resources/coloc/{study}",
        dir_bfile  = "resources/genotypes/{study}",
        gene_loc = ancient("resources/misc/gencode.v44.gene_type.tsv")
    output:
        cmd_file = "cmds/mvcoloc/{study}/main_cmds.txt",
        rerun_file = "cmds/mvcoloc/{study}/rerun_cmds.txt"
    params:
        script = "workflow/rules/snakescripts/coloc/run_mvcoloc_cli.R",
        traits = lambda wc: " ".join(TRAITS),
        outdir = "results/mvcoloc/{study}",
        window_bp = 100000,
        runsusie_coverage = 0.1,
        min_p_gwas = 1e-4,
        runsusie_maxit = 200,
        runsusie_repeat = "FALSE",
        runsusie_timeout = 180,
        coloc_timeout = 180,
        p12 = 1e-5,
        nthreads = 1
    conda: "../envs/renv.yaml"
    shell:
        r"""
        PROJECT_ROOT=$(pwd)
        CONDA_ENV="$CONDA_PREFIX"
        SCRIPT_DIR=$(dirname {output.cmd_file})

        mkdir -p "$SCRIPT_DIR" > {output.cmd_file}
        mkdir -p "$SCRIPT_DIR" > {output.rerun_file}

        for chr in $(seq 1 22); do
            for pheno in {params.traits}; do
                outfile="{params.outdir}/${{pheno}}/chr${{chr}}.mvcoloc.tsv"
                log="{params.outdir}/${{pheno}}/chr${{chr}}.mvcoloc.tsv.log"
                script_file="$SCRIPT_DIR/task_chr${{chr}}_${{pheno}}.sh"
                susie_rds="results/susie_gwas/{wildcards.study}/${{pheno}}/chr${{chr}}.susie.rds"
                cat > "$script_file" <<TASK_EOF
#!/bin/bash
set -euo pipefail
cd $PROJECT_ROOT
if [ -f "$outfile" ] && [ -s "$outfile" ]; then
    echo "\$(date): SKIP $outfile"
    exit 0
fi

export PATH=$CONDA_ENV/bin:\$PATH
export CONDA_PREFIX=$CONDA_ENV
export R_LIBS_SITE=$CONDA_ENV/lib/R/library

mkdir -p \$(dirname $outfile)

Rscript {params.script} \
    --study {wildcards.study} \
    --chr $chr \
    --pheno $pheno \
    --ld_eqtl {input.dir_ld_eqtl}/chr${{chr}} \
    --dir_eqtl {input.dir_eqtl} \
    --susie_gwas $susie_rds \
    --gene_loc {input.gene_loc} \
    --dir_bfile {input.dir_bfile} \
    --output $outfile \
    --threads {params.nthreads} \
    --window_bp {params.window_bp} \
    --p12 {params.p12} \
    --runsusie_coverage {params.runsusie_coverage} \
    --runsusie_maxit {params.runsusie_maxit} \
    --runsusie_repeat {params.runsusie_repeat} \
    --runsusie_timeout {params.runsusie_timeout} \
    --coloc_timeout {params.coloc_timeout} \
    > $log 2>&1
TASK_EOF

                chmod +x "$script_file"
                if [ ! -f "$outfile" ]; then
                    echo "$PROJECT_ROOT/$script_file" >> {output.rerun_file}
                else
                    echo "$PROJECT_ROOT/$script_file" >> {output.cmd_file}
                fi
            done
        done

        echo "Generated $(wc -l < {output.cmd_file}) main cmd scripts"
        echo "Generated $(wc -l < {output.rerun_file}) rerun cmd scripts"
        """


rule run_mvcoloc:
    """
    Execute multivariant coloc tasks via nci-parallel.
    """
    input:
        cmd_file = "cmds/mvcoloc/{study}/{batch}_cmds.txt"
    output:
        done = touch("results/mvcoloc/{study}/.done.{batch}")
    resources:
        ncpus  = 96,
        mem    = "380G",
        time   = "24:00:00",
        jobfs  = "4G"
    params:
        ncores_per_task        = 1,  
        ncores_per_numanode     = 12, # 48  → ppr
        timeout         = 10000
    log: "logs/mvcoloc/{study}/nciparallel.{batch}.log"
    shell:
        """
        module load nci-parallel/1.0.0a 2>> {log}

        NCMDS=$(wc -l < {input.cmd_file})

        M=$(( PBS_NCI_NCPUS_PER_NODE / {params.ncores_per_task} )) 
        PPR=$(( {params.ncores_per_numanode} / {params.ncores_per_task} ))
        NP=$((M * PBS_NCPUS / PBS_NCI_NCPUS_PER_NODE))
        
        echo $PPR
        echo $NP
        mpirun -np $NP \
            --map-by ppr:${{PPR}}:NUMA:PE={params.ncores_per_task} \
            nci-parallel \
            --verbose \
            --input-file {input.cmd_file} \
        >> {log} 2>&1 \
        && touch {output.done}
        """

rule run_mv_coloc_debug:
    """
    Debug: run a single chr × trait on 1 node without nci-parallel.
    Must be invoked explicitly:
        snakemake --allowed-rules run_mv_coloc_debug results/mvcoloc/.../debug.mvcoloc.tsv
    """
    input:
        dir_eqtl = ancient("resources/coloc/{study}"),
        ld_eqtl = ancient("resources/ld/eqtl/{study}/chr{chr}"),
        gwas = "resources/ma/{pheno}.ma",
        gene_loc = ancient("resources/misc/gencode.v44.gene_type.tsv"),
        pheno_metadata = ancient("metadata/trait.tsv"),
        dir_bfile = ancient("resources/genotypes/{study}")
    output:
        coloc = "results/mvcoloc/{study}/chr{chr}/{pheno}/debug.mvcoloc.tsv"
    threads: 1
    log: "logs/mvcoloc/{study}/chr{chr}/{pheno}/mvcoloc_debug.log"
    resources:
        mem = "16G",
        ncpus = 1,
        jobfs = "4G",
        time = "24:00:00"
    params:
        window_bp = 100000,
        runsusie_coverage = 0.1,
        p12 = 1e-5,
        runsusie_maxit = 200,
        runsusie_repeat = False
    conda: "../envs/renv.yaml"
    script: "snakescripts/coloc/run_multivariant_coloc.R"


def _target_mvcoloc_tsv(x):
    with open("resources/misc/target_phenotypes.txt") as f:
        PHENOS = [line.strip() for line in f if line.strip()]
    return expand(
        f"results/mvcoloc/{x.study}/{{pheno}}/chr{{chr}}.mvcoloc.tsv",
        chr=range(1, 23), pheno=PHENOS
    )


rule concat_mv_coloc_all:
    """Concatenate all mvcoloc results into a single parquet file."""
    input: _target_mvcoloc_tsv
    output:
        "results/aggregate/coloc/{study}.mvcoloc.parquet.gz"
    conda: "../envs/renv.yaml"
    log: "logs/aggregate/{study}.concat_mvcoloc.log"
    threads: 8
    resources:
        mem = "64G",
        ncpus = 8
    script: "snakescripts/aggregate/concat_mvcoloc_all.R"
