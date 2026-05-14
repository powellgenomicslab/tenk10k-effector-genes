# Rules to run scDRS

wildcard_constraints:
    # exclude forward slash
    study="[^/]+",
    phenotype="[^/]+"

rule magma_to_zscore:
    """
    Create zscore file for scDRS from collection of magma file
    format: https://martinjzhang.github.io/scDRS/file_format.html#pval-file-zscore-file
    """
    input: "results/aggregate/{study}.magma.gz.parquet"
    output: "resources/scdrs/zscore/{study}.zscore.tsv"
    conda: "renv"
    log: "logs/scdrs/{study}.magma_to_zscore.log"
    script: "snakescripts/scdrs/magma_to_zscore.R"

rule scdrs_munge_gs:
    """
    Create .gs file from .zscore file for scDRS
    """
    input: "resources/scdrs/zscore/{study}.zscore.tsv"
    output: "resources/scdrs/gs/{study}.gs"
    conda: "scverse"
    log: "logs/scdrs/{study}.scdrs_munge_gs.log"
    params:
        fdr = 0.05,
        n_max = 1000
    shell:
        """
        scdrs munge-gs \
            --out-file {output} \
            --zscore-file {input} \
            --weight zscore \
            --n-max {params.n_max} \
            --fdr {params.fdr}
        """

# def get_covar_study(x):
#     dir_covar = Path(f"resources/scdrs/sc_covar/{x.study}")
#     return [str(f) for f in list(dir_covar.glob("*.covar.csv"))]

rule h5ad_annot:
    """
    Convert .h5ad file to .rds file for scPred annotation
    """
    input:
        h5ad = "resources/scdrs/h5ad/{study}.prep.h5ad",
        ref_scpred = "resources/scpred/{study}.{annot}.ref.rds"
    output:
        raw = "results/scdrs/annot_sample/{study}.{annot}.raw.csv",
        format = "results/scdrs/annot_sample/{study}.{annot}.format.csv"
    conda: "r-sc"
    resources:
        ncpus = 8,
        mem =  "300G",
        queue = "hugemem",
        jobfs = "20G"
    log: "logs/scdrs/{study}.{annot}.h5ad_to_rds.log"
    script: "snakescripts/scdrs/annot/{wildcards.annot}.R"

checkpoint chunk_gs:
    """
    Chunk .gs file into smaller files (per phenotype) for scDRS
    """
    input:
        gs = "resources/scdrs/gs/{study}.gs",
        config = "resources/scdrs/config/{study}.yaml"
    output:
        directory("resources/scdrs/gs_chunked/{study}")
    conda: "scverse"
    log: "logs/scdrs/{study}.chunk_gs.log"
    shell:
        """
        mkdir -p {output}
        awk -v DIR={output} \
            'NR == 1 {{HEAD = $0; next}}
             {{TRAIT = $1;
              print HEAD RS $0 > DIR "/" TRAIT ".gs"}}' \
            {input.gs}
        """

rule scdrs_prep_h5ad_cov:
    """
    Prepare .h5ad file for scDRS
    """
    input:
        h5ad = "resources/scdrs/h5ad/{study}.h5ad",
        sample_covar = "resources/scdrs/sample_covar/{study}.sample_covar.csv",
        config = "resources/scdrs/config/{study}.yaml"
    output:
        prep_h5ad = "resources/scdrs/h5ad/{study}.prep.h5ad",
        cov = "resources/scdrs/cov/{study}.cov.tsv"
    conda: "scverse"
    resources:
        ncpus = 16,
        mem =  "480G",
        queue = "hugemem"
    log: "logs/scdrs/{study}.prep_h5ad_cov.log"
    script: "snakescripts/scdrs/prep/{wildcards.study}.py"


rule scdrs_regress_h5ad_cov:
    """
    Prepare regressed .h5ad file for scDRS
    """
    input:
        prep_h5ad = "resources/scdrs/h5ad/{study}.prep.h5ad",
        # gs = "resources/scdrs/gs/{study}.gs",
        cov = "resources/scdrs/cov/{study}.cov.tsv",
        config = "resources/scdrs/config/{study}.yaml"
    output:
        # save as pickle as cannot save with standard .h5ad
        "resources/scdrs/h5ad/{study}.reg.h5ad.pkl"
    conda: "scverse"
    resources:
        ncpus = 16,
        mem =  "64G"
    log: "logs/scdrs/{study}.regress_h5ad_cov.log"
    script: "snakescripts/scdrs/regress_h5ad_cov.py"

rule scdrs_compute_score_api:
    """
    Compute scDRS score using the python API per phenotype
    """
    input:
        reg_pkl = "resources/scdrs/h5ad/{study}.reg.h5ad.pkl",
        gs = "resources/scdrs/gs_chunked/{study}/{phenotype}.gs",
        config = "resources/scdrs/config/{study}.yaml"
    output:
        cell = "results/scdrs/cell_score/{study}/{phenotype}.cell_score.tsv.parquet.gz",
        celltype = "results/scdrs/cell_type_stats/{study}/{phenotype}.cell_type_stats.tsv"
    conda: "scverse"
    resources:
        ncpus = 8,
        mem =  "64G"
    log: "logs/scdrs/{study}-{phenotype}.regress_h5ad_cov.log"
    script: "snakescripts/scdrs/compute_score.py"

rule scdrs_get_top_score:
    """
    Calculate counts of cells annotate with specific cell types
    within the top nth quantile based on the normalised scDRS score
    """
    input:
        reg_pkl = "resources/scdrs/h5ad/{study}.reg.h5ad.pkl",
        cell_score = "results/scdrs/cell_score/{study}/{phenotype}.cell_score.tsv.parquet.gz",
        config = "resources/scdrs/config/{study}.yaml"
    output: "results/scdrs/cell_type_top/{study}/{phenotype}.cell_type_top.tsv"
    conda: "scverse"
    params:
        percentiles = [90, 95,99]
    resources:
        ncpus = 8,
        mem =  "48G"
    log: "logs/scdrs/{study}-{phenotype}.get_top_score.log"
    script: "snakescripts/scdrs/get_top_score.py"

rule scdrs_annot_stats:
    """
    Compute cell-type level stats for custom annotations
    """
    input:
        reg_pkl = "resources/scdrs/h5ad/{study}.reg.h5ad.pkl",
        cell_score = "results/scdrs/cell_score/{study}/{phenotype}.cell_score.tsv.parquet.gz",
        config = "resources/scdrs/config/{study}.yaml",
        annot = "results/scdrs/annot_sample/{study}.{annot}.format.csv"
    output: "results/scdrs/cell_type_annot/{study}/{phenotype}.{annot}.stats.tsv"
    conda: "scverse"
    resources:
        ncpus = 8,
        mem =  "48G"
    log: "logs/scdrs/{study}-{phenotype}-{annot}.annot_stats.log"
    script: "snakescripts/scdrs/cell_type_stats_annot.py"



# aggregate scdrs per phenotype
def scdrs_aggregate_stat(x):
    """
    Aggregate scDRS cell type stats per phenotype
    """
    dir_gs = Path(str(checkpoints.chunk_gs.get(**x).output[0]))
    gs_files = dir_gs.glob("*.gs")
    traits = [str(gs_file.stem) for gs_file in gs_files]

    celltype_stats_files = [f"results/scdrs/cell_type_stats/{x.study}/{p}.cell_type_stats.tsv" for p in traits]
    celltype_top_files = [f"results/scdrs/cell_type_top/{x.study}/{p}.cell_type_top.tsv" for p in traits]
    return {'cell_type_stats':celltype_stats_files,
            'cell_type_top': celltype_top_files}

rule scdrs_aggregate_stat:
    input: unpack(scdrs_aggregate_stat)
    output:
        cell_type_stats = "results/aggregate/{study}.scdrs.cell_type_stats.tsv",
        cell_type_top = "results/aggregate/{study}.scdrs.cell_type_top.tsv"
    resources:
        ncpus = 8,
        mem =  "64G"
    log: "logs/aggregate/{study}.scdrs_cell_type_stat.log"
    conda: "renv"
    script: "snakescripts/aggregate/scdrs_cell_type.R"
    

def scdrs_aggregate_stat_annot(x):
    """
    Aggregate scDRS cell type stats per phenotype for custom annotations
    """
    dir_gs = Path(str(checkpoints.chunk_gs.get(**x).output[0]))
    gs_files = dir_gs.glob("*.gs")
    traits = [str(gs_file.stem) for gs_file in gs_files]

    celltype_stats_files = [f"results/scdrs/cell_type_annot/{x.study}/{p}.{x.annot}.stats.tsv" for p in traits]
    return celltype_stats_files

rule scdrs_aggregate_stat_annot:
    input: scdrs_aggregate_stat_annot
    output: "results/aggregate/{study}.scdrs.{annot}.stats.tsv",
    resources:
        ncpus = 8,
        mem =  "64G"
    log: "logs/aggregate/{study}.scdrs.{annot}.stats.log"
    conda: "renv"
    script: "snakescripts/aggregate/scdrs_custom_annot.R"

def scdrs_aggregate_score(x):
    """
    Aggregate scDRS cell-level scores per phenotype
    """
    dir_gs = Path(str(checkpoints.chunk_gs.get(**x).output[0]))
    gs_files = dir_gs.glob("*.gs")
    traits = [str(gs_file.stem) for gs_file in gs_files]

    cell_score_files = [f"results/scdrs/cell_score/{x.study}/{p}.cell_score.tsv.parquet.gz" for p in traits]
    return cell_score_files

rule scdrs_aggregate_score:
    """
    Get umap coordinates (and annotations) of sampled cells
    and cell-level scores for plotting
    """
    input:
        cell_score = scdrs_aggregate_score,
        prep_h5ad = "resources/scdrs/h5ad/{study}.prep.h5ad",
        config = "resources/scdrs/config/{study}.yaml"
    output:
        scores  = "results/aggregate/{study}.scdrs.cell_score.tsv.parquet.gz",
        mcp  = "results/aggregate/{study}.scdrs.cell_mcp.tsv.parquet.gz"
    conda: "scverse"
    resources:
        ncpus = 16,
        mem =  "80G"
    log: "logs/aggregate/{study}.scdrs_cell_score.log"
    script: "snakescripts/aggregate/scdrs_cell_score.py"



# rule scdrs_prep_cov:
#     """
#     Prepare covariates file from .h5ad file
#     """
#     input:
#         h5ad = "resources/scdrs/h5ad/{study}.prep.h5ad",
#         # gs = "resources/scdrs/gs/{geno_set}.gs",
#         # config = "resources/scdrs/config/{study}.yaml",
#         sample_covar = "resources/scdrs/sample_covar/{study}.sample_covar.csv"
#     output:
#         cov = "resources/scdrs/cov/{study}.cov.tsv"
#         # regressed_h5ad = "resources/scdrs/h5ad/{geno_set}-{study}.regressed.h5ad"
#     conda: "scverse"
#     log: "logs/scdrs/prep_cov.{study}.log"
#     resources:
#         ncpus = 16,
#         mem =  "600G",
#         time = "24:00:00",
#         queue = "hugemem"
#     script: "snakescripts/scdrs/prep_cov.py"


# def scdrs_read_params(x, key):
#     with open(f"resources/scdrs/config/{x.study}.yaml", 'r') as f:
#         config = yaml.full_load(f)
#     return config[key]


# rule scdrs_compute_score:
#     """
#     Extract covariate (.cov) file from .h5ad file for scDRS
#     """
#     input:
#         gs = "resources/scdrs/gs/{study}.gs",
#         prep_h5ad = "resources/scdrs/h5ad/{study}.prep.h5ad",
#         cov = "resources/scdrs/cov/{study}.cov.tsv",
#         config = "resources/scdrs/config/{study}.yaml"
#     output: directory("results/scdrs/{study}")
#     conda: "scverse"
#     resources:
#         ncpus = 24,
#         mem =  "600G",
#         queue = "hugemem",
#         time = "48:00:00"
#     log: "logs/scdrs/compute_score.{study}.log"
#     params:
#         h5ad_species = lambda x: scdrs_read_params(x, 'h5ad_species'),
#         gs_species = lambda x: scdrs_read_params(x, 'gs_dst_species'),
#         filter_data = lambda x: scdrs_read_params(x, 'filter_data'),
#         raw_count = lambda x: scdrs_read_params(x, 'raw_count'),
#         n_ctrl = lambda x: scdrs_read_params(x, 'n_ctrl'),
#         adj_prop_col = lambda x: scdrs_read_params(x, 'adj_prop_col'),
#         return_ctrl_raw_score = lambda x: scdrs_read_params(x, 'return_ctrl_raw_score'),
#         return_ctrl_norm_score = lambda x: scdrs_read_params(x, 'return_ctrl_norm_score')
#     shell:
#         """
#         mkdir -p {output}
#         scdrs compute-score \
#             --h5ad-file {input.prep_h5ad} \
#             --h5ad-species {params.h5ad_species} \
#             --gs-file {input.gs} \
#             --gs-species {params.gs_species} \
#             --out-folder {output} \
#             --cov-file {input.cov} \
#             --flag-filter-data {params.filter_data} \
#             --flag-raw-count {params.raw_count} \
#             --n-ctrl {params.n_ctrl} \
#             --adj-prop {params.adj_prop_col} \
#             --flag-return-ctrl-raw-score {params.return_ctrl_raw_score}  \
#             --flag-return-ctrl-norm-score {params.return_ctrl_norm_score}  
#         """

