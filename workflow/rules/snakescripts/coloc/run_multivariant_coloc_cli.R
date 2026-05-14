# CLI wrapper for nci-parallel. Builds a mock snakemake list and sources
# the main run_multivariant_coloc.R script.
#
# Usage: Rscript run_multivariant_coloc_cli.R --study <s> --chr <c> --pheno <p> ...

# ------------------------------------------------------------------
# Lock down threading BEFORE anything else — must be first!
# ------------------------------------------------------------------
Sys.setenv(
    OMP_NUM_THREADS         = "1",
    MKL_NUM_THREADS         = "1",
    OPENBLAS_NUM_THREADS    = "1",
    VECLIB_MAXIMUM_THREADS  = "1",
    R_DATATABLE_NUM_THREADS = "1"
)

library(optparse)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--study",             type = "character"),
    make_option("--chr",               type = "character"),
    make_option("--pheno",             type = "character"),
    make_option("--dir_eqtl",          type = "character"),
    make_option("--ld_eqtl",           type = "character"),
    make_option("--gwas",              type = "character"),
    make_option("--gene_loc",          type = "character"),
    make_option("--pheno_metadata",    type = "character"),
    make_option("--dir_bfile",         type = "character"),
    make_option("--output",            type = "character"),
    make_option("--threads",           type = "integer",  default = 48),
    make_option("--window_bp",         type = "integer",  default = 100000),
    make_option("--runsusie_coverage", type = "double",   default = 0.1),
    make_option("--p12",               type = "double",   default = 1e-5),
    make_option("--runsusie_maxit",    type = "integer",  default = 200),
    make_option("--runsusie_repeat",   type = "logical",  default = FALSE),
    make_option("--mem_limit_gb",      type = "integer",  default = 16)
)))

# Validate required args
required <- c("study", "chr", "pheno", "dir_eqtl", "ld_eqtl", "gwas",
              "gene_loc", "pheno_metadata", "dir_bfile", "output")
missing <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

cat(sprintf("[%s] PID %d | %s chr%s %s | threads=%d\n",
            Sys.time(), Sys.getpid(), opt$study, opt$chr, opt$pheno, opt$threads))

# Print all resolved paths for debugging
cat(sprintf("  dir_eqtl:      %s (exists=%s)\n", opt$dir_eqtl,      dir.exists(opt$dir_eqtl)))
cat(sprintf("  ld_eqtl:       %s (exists=%s)\n", opt$ld_eqtl,       dir.exists(opt$ld_eqtl)))
cat(sprintf("  gwas:          %s (exists=%s)\n", opt$gwas,           file.exists(opt$gwas)))
cat(sprintf("  gene_loc:      %s (exists=%s)\n", opt$gene_loc,       file.exists(opt$gene_loc)))
cat(sprintf("  pheno_metadata:%s (exists=%s)\n", opt$pheno_metadata, file.exists(opt$pheno_metadata)))
cat(sprintf("  dir_bfile:     %s (exists=%s)\n", opt$dir_bfile,      dir.exists(opt$dir_bfile)))
cat(sprintf("  output:        %s\n", opt$output))
cat(sprintf("  threads:       %d\n", opt$threads))

# ------------------------------------------------------------------
# Build mock snakemake object
# ------------------------------------------------------------------
snakemake <- list(
    input = list(
        dir_eqtl       = opt$dir_eqtl,
        ld_eqtl        = opt$ld_eqtl,
        gwas           = opt$gwas,
        gene_loc       = opt$gene_loc,
        pheno_metadata = opt$pheno_metadata,
        dir_bfile      = opt$dir_bfile
    ),
    output = list(
        coloc = opt$output
    ),
    wildcards = list(
        study = opt$study,
        pheno = opt$pheno,
        chr   = opt$chr
    ),
    params = list(
        window_bp         = opt$window_bp,
        runsusie_coverage = opt$runsusie_coverage,
        p12               = opt$p12,
        runsusie_maxit    = opt$runsusie_maxit,
        runsusie_repeat   = opt$runsusie_repeat,
        mem_limit_gb     = opt$mem_limit_gb
    ),
    threads = opt$threads
)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Source the main script
# ------------------------------------------------------------------
tryCatch({
    script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
    main_script <- if (!is.null(script_dir)) {
        file.path(script_dir, "run_multivariant_coloc.R")
    } else {
        "workflow/rules/snakescripts/coloc/run_multivariant_coloc.R"
    }

    if (!file.exists(main_script)) {
        stop("Main script not found: ", main_script)
    }

    cat(sprintf("  Sourcing: %s\n", main_script))
    source(main_script)
}, error = function(e) {
    cat(sprintf("[%s] FAILED: %s | %s\n", Sys.time(), opt$output, e$message))
    traceback()
    quit(status = 1)
})