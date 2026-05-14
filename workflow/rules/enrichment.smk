# Perform gene set enrichment analysis of SMR results
# use gget enrichr api
# https://pachterlab.github.io/gget/en/enrichr.html

rule gget_enrichr:
    """Enrichr enrichment analysis using gget
    note: this requires internet access to connect to external databases
    """
    input:
        gene_universe = "results/enrichment/{study}/gene_universe.txt",
        msmr_sig = "results/aggregate/msmr_sig/{study}~q_{q_thresh}~heidi_{heidi_thresh}.msmr_sig.tsv",
        # full gene_set list: https://maayanlab.cloud/Enrichr/#libraries
        gene_set = "resources/misc/enrichr.gene_set.txt"
    conda:
        "pydata"
    output:
        phenotype = "results/enrichment/{study}~q_{q_thresh}~heidi_{heidi_thresh}.enrichr.phenotype.tsv",
        biosample = "results/enrichment/{study}~q_{q_thresh}~heidi_{heidi_thresh}.enrichr.biosample.tsv"
    script: "snakescripts/enrichment/gget_enrichr.py"

checkpoint gget_enrichr_dir_pheno:
    """Enrichr enrichment analysis per phenotype using gget
    note: this requires internet access to connect to external databases
    """
    input: msmr_sig = "results/aggregate/msmr_sig/{study}~q_{q_thresh}~heidi_{heidi_thresh}.msmr_sig.tsv"
    log: "logs/enrichment/gget_enrichr_dir_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}.log"
    output:
        dir_pheno = directory("results/enrichment_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/gene_celltype")
    shell:
        """
        mkdir -p {output.dir_pheno}
        cut -f1,2,5 {input.msmr_sig} | \
            sort -u | \
            awk 'NR > 1 {{print $1,$3 > "{output.dir_pheno}/" $2 ".txt"}}'
        """

rule gget_enrichr_pheno:
    """Enrichr enrichment analysis per phenotype using gget
    note: this requires internet access to connect to external databases
    """
    input:
        gene_universe = "results/enrichment/{study}/gene_universe.txt",
        msmr_sig = "results/enrichment_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/gene_celltype/{pheno}.txt",
        # full gene_set list: https://maayanlab.cloud/Enrichr/#libraries
        gene_set = "resources/misc/enrichr.gene_set.txt"
    conda: "pydata"
    log: "logs/enrichment/gget_enrichr_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/{pheno}.log"
    params: min_gene = 5
    resources:
        queue = "copyq",
        ncpus = "1",
        mem = "4GB",
        time = "06:00:00"
    output:
        pheno = "results/enrichment_pheno/{study}~q_{q_thresh}~heidi_{heidi_thresh}/enrichr/{pheno}.tsv"
    script: "snakescripts/enrichment/gget_enrichr_pheno.py"

def gget_enrichr_pheno_aggregate(x):
    DIR = Path(checkpoints.gget_enrichr_dir_pheno.get(**x).output[0])
    pheno = [f.with_suffix('').name for f in DIR.glob("*.txt")]
    return [f"results/enrichment_pheno/{x.study}~q_{x.q_thresh}~heidi_{x.heidi_thresh}/enrichr/{p}.tsv"
            for p in pheno]

rule aggregate_gget_enrichr:
    """Aggregate Enrichr enrichment analysis per phenotype using gget
    note: this requires internet access to connect to external databases
    """
    input: gget_enrichr_pheno_aggregate
    output: "results/aggregate/{study}/enrichr.q_{q_thresh}~heidi_{heidi_thresh}.tsv.gz"
    log: "logs/aggregate/enrichr.{study}.q_{q_thresh}.heidi_{heidi_thresh}.log"
    shell: "awk 'NR == FNR || FNR > 1' {input} | gzip -c > {output}"

# new enrichment analysis with gprofiler
checkpoint prep_gene_set:
    wildcard_constraints:
        set = "[^/]+"
    output: directory("resources/enrichment/set/{set}")
    conda: "renv"
    log: "logs/prep_gene_set/{set}.log"
    resources:
        ncpus = 8,
        mem = "64G",
        time = "06:00:00"
    script: "snakescripts/prep_gene_set/{wildcards.set}.R"

rule map_string_ids:
    """Map gencode Ensembl gene IDs / HGNC symbols to STRING protein IDs.
    One-time resource generation; output is reused by all prep_gene_set_string scripts.
    note: this requires internet access to query the STRING API
    """
    input:
        gencode = "resources/misc/gencode.v44.gene_type.tsv"
    conda: "pydata"
    output:
        mapping = "resources/enrichment/string_id_map.tsv"
    params:
        species = 9606,
        string_version = "12.0",
        caller_identity = "tenk10k_smr",
        batch_size = 2000,
        retries = 3
    resources:
        ncpus = "1",
        queue = "copyq",
        mem = "4GB",
        time = "02:00:00"
    log: "logs/enrichment/map_string_ids.log"
    script: "snakescripts/enrichment/map_string_ids.py"

checkpoint prep_gene_set_string:
    wildcard_constraints:
        set = "[^/]+"
    input:
        string_id_map = "resources/enrichment/string_id_map.tsv"
    output: directory("resources/enrichment/set_string/{set}")
    conda: "renv"
    log: "logs/prep_gene_set_string/{set}.log"
    resources:
        ncpus = 8,
        mem = "64G",
        time = "06:00:00"
    script: "snakescripts/prep_gene_set_string/{wildcards.set}.R"

# --- Enrichment methods dispatch ---
# Explicit method lists — add new methods to the right list
PYTHON_ENRICH_METHODS = ["string_rank_api"]
R_ENRICH_METHODS      = ["gprofiler", "stringdb"]

ruleorder: enrichment_by_method_py > enrichment_by_method

rule enrichment_by_method_py:
    """Enrichment analysis using Python-based methods
    (e.g. STRING values/ranks enrichment API)
    note: this requires internet access to connect to external databases
    """
    wildcard_constraints:
        enrich_method = "|".join(PYTHON_ENRICH_METHODS)
    input:
        dir_gene_set = "resources/enrichment/set_string/{set}/{biosample}",
        gene_universe = "results/enrichment/{study}/gene_universe.txt"
    conda: "pydata"
    output:
        enrich    = "results/enrichment/{study}/{set}/{biosample}.{enrich_method}.tsv",
        ppi_network = "results/enrichment/{study}/{set}/{biosample}.{enrich_method}.ppi_network.tsv",
        bipartite   = "results/enrichment/{study}/{set}/{biosample}.{enrich_method}.bipartite.tsv"
    params:
        # --- shared params ---
        min_genes = 5,
        # --- string values/ranks API params ---
        species = 9606,
        string_version = "12.0",
        fdr_threshold = 0.05,
        caller_identity = "tenk10k_smr",
        poll_interval = 30,      # seconds between status checks
        max_wait = 7200,         # max total wait time (2 h)
        api_key_file = ".string_api_key",
        premap_ids = False,      # set True to pre-map gene symbols → STRING IDs
        score_threshold = 400,   # min combined score for PPI network edges
        rank_direction = 0       # 0=both extremes, +1=top only, -1=bottom only
    resources:
        ncpus = "1",
        queue = "copyq",
        mem = "4GB",
        time = "10:00:00"
    log: "logs/enrichment/{enrich_method}/{study}/{set}/{biosample}.{enrich_method}.log"
    script: "snakescripts/enrichment/{wildcards.enrich_method}.py"

rule enrichment_by_method:
    """R-based enrichment analysis (gprofiler, stringdb, etc.)
    note: this requires internet access to connect to external databases
    """
    wildcard_constraints:
        enrich_method = "|".join(R_ENRICH_METHODS)
    input:
        # gene_set = "resources/enrichment/set/{set}/{biosample}/{pheno}.txt",
        dir_gene_set = "resources/enrichment/set/{set}/{biosample}",
        gene_universe = "results/enrichment/{study}/gene_universe.txt"
    conda: "renv"
    output: enrich = "results/enrichment/{study}/{set}/{biosample}.{enrich_method}.tsv"
    params:
        # --- gprofiler params ---
        sources = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", "HP", "HPA", "TRANSFAC", "CORUM"],
        # GO BP 2024 gmt file
        custom_gmt = {"GO_BP_2024": "gp__vyDj_re4A_Cx4"},
        domain_scope = "custom",
        highlight = True,
        ordered_query = False,
        multi_query = False,
        significant = True,
        # --- shared params ---
        min_genes = 5,
        # --- stringdb params ---
        species = 9606,
        score_threshold = 400,
        string_version = "12.0",
        fdr_threshold = 0.05
    resources:
        # opt = "-l ood=jupyter ",
        ncpus = "1",
        queue = "copyq",
        mem = "4GB",
        time = "06:00:00"
    log: "logs/enrichment/{enrich_method}/{study}/{set}/{biosample}.{enrich_method}.log"
    script: "snakescripts/enrichment/{wildcards.enrich_method}.R"

def _get_enrich(x):
    if x.enrich_method in PYTHON_ENRICH_METHODS:
        DIR = checkpoints.prep_gene_set_string.get(**x).output[0]
    else:
        DIR = checkpoints.prep_gene_set.get(**x).output[0]
    BIOSAMPLES = [f.name for f in Path(DIR).glob("*") if f.is_dir()]
    enrich_files = [f"results/enrichment/{x.study}/{x.set}/{b}.{x.enrich_method}.tsv" for b in BIOSAMPLES]
    return enrich_files

rule aggregate_enrichment:
    """Aggregate enrichment results
    """
    input: _get_enrich
    output: "results/aggregate/enrichment/{study}.{set}.{enrich_method}.tsv.gz"
    log: "logs/aggregate/enrichment/{study}.{set}.{enrich_method}.log"
    shell:
        """
        (
            awk -v OFS='\\t' 'FNR == 1 {{print "cell_type", $0; next}}' $(echo {input} | tr ' ' '\n' | head -1)
            for f in {input}; do
                cell=$(basename "$f" | sed 's/\\..*$//')
                awk -v cell="$cell" -v OFS='\\t' 'FNR > 1 {{print cell, $0}}' "$f"
            done
        ) | gzip -c > {output}
        """
