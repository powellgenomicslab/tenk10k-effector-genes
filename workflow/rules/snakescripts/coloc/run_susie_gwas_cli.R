# CLI wrapper for nci-parallel
library(optparse)
library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(sys)
library(susieR)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--study",             type = "character"),
    make_option("--chr",               type = "character"),
    make_option("--pheno",             type = "character"),
    make_option("--gwas",              type = "character"),
    make_option("--gene_loc",          type = "character"),
    make_option("--pheno_metadata",    type = "character"),
    make_option("--dir_bfile",         type = "character"),
    make_option("--output",            type = "character"),
    make_option("--threads",           type = "integer",  default = 8),
    make_option("--window_bp",         type = "integer",  default = 100000),
    make_option("--runsusie_coverage", type = "double",   default = 0.1),
    make_option("--min_p_gwas",        type = "double",   default = 1e-4),
    make_option("--runsusie_maxit",    type = "integer",  default = 200),
    make_option("--runsusie_repeat",   type = "logical",  default = FALSE))
))

# Validate required args
required <- c("study", "chr", "pheno", "gwas", "gene_loc", "pheno_metadata", "dir_bfile", "output")
missing <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

cat(sprintf("[%s] PID %d | %s chr%s %s | threads=%d\n",
            Sys.time(), Sys.getpid(), opt$study, opt$chr, opt$pheno, opt$threads))

# Print all resolved paths for debugging
cat(sprintf("  gwas:          %s (exists=%s)\n", opt$gwas,           file.exists(opt$gwas)))
cat(sprintf("  gene_loc:      %s (exists=%s)\n", opt$gene_loc,       file.exists(opt$gene_loc)))
cat(sprintf("  pheno_metadata:%s (exists=%s)\n", opt$pheno_metadata, file.exists(opt$pheno_metadata)))
cat(sprintf("  dir_bfile:     %s (exists=%s)\n", opt$dir_bfile,      dir.exists(opt$dir_bfile)))
cat(sprintf("  output:        %s\n", opt$output))
cat(sprintf("  threads:       %d\n", opt$threads))

# ------------------------------------------------------------------
# Assign options
# ------------------------------------------------------------------

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)


# ==========================================================================
# 1. Parameters
# ==========================================================================
gwas_file         <- opt$gwas
gene_loc_file     <- opt$gene_loc
pheno_metadata_file <- opt$pheno_metadata
dir_bfile         <- opt$dir_bfile
output_file       <- opt$output

window_bp <- opt$window_bp
n_threads <- opt$threads
min_p_gwas <- opt$min_p_gwas

pheno <- opt$pheno
study <- opt$study
chr   <- opt$chr

# Progress log (visible from outside future workers)
progress_log <- paste0(output_file, ".progress.log")
file.create(progress_log)

log_progress <- function(msg) {
    tryCatch(
        write(paste0("[", Sys.time(), "] [PID:", Sys.getpid(), "] ", msg),
              file = progress_log, append = TRUE),
        error = function(e) invisible(NULL)
    )
}

log_progress(glue("Start: study={study} chr={chr} pheno={pheno} threads={n_threads}"))

io_threads <- min(n_threads, 8)
setDTthreads(io_threads)

cat(glue("Loading data with {io_threads} threads..."), "\n")

# --- Phenotype metadata ---
t0 <- Sys.time()
df_pheno_meta <- fread(pheno_metadata_file)[include == TRUE & trait_id == pheno]
pheno_type <- ifelse(df_pheno_meta$supercategory == "biological", "quant", "cc")
pheno_sd   <- ifelse(pheno_type == "quant", 1, NA)
cat(glue("  Pheno metadata: {round(difftime(Sys.time(), t0, units='secs'), 1)}s"), "\n")

# --- Gene locations ---
t0 <- Sys.time()
CHR <- chr
df_gene_loc <- fread(gene_loc_file)[chr == CHR]
df_gene_loc[, `:=`(cis_start = pmax(1L, start - window_bp), cis_end = end + window_bp)]
cat(glue("  Gene locations: {nrow(df_gene_loc)} genes ({round(difftime(Sys.time(), t0, units='secs'), 1)}s)"), "\n")

# --- GWAS — load full chr, keep SNPs seen in any eQTL file ---
t0 <- Sys.time()
df_gwas <- fread(gwas_file)
dup_snps <- df_gwas[duplicated(SNP), unique(SNP)]
df_gwas <- df_gwas[!SNP %in% dup_snps]
df_gwas[, position := str_split(SNP, ":") |> map_chr(2) |> as.integer()]

cat(glue("  GWAS: {nrow(df_gwas)} rows ({round(difftime(Sys.time(), t0, units='secs'), 1)}s)"), "\n")

# filter only region with egene and P_gwas < threshold
df_gwas_top_pos <- df_gwas[P < min_p_gwas, .(position)]
matched_indices <- df_gene_loc[df_gwas_top_pos, on = .(cis_start <= position, cis_end >= position), which = TRUE, nomatch = 0]
gwas_genes <- unique(df_gene_loc[matched_indices])$ensembl_gene_id

# ==========================================================================
# 3. Run Coloc SuSiE
# ==========================================================================
jobfs_dir <- Sys.getenv("PBS_JOBFS", tempdir())
cat(glue("  jobfs_dir: {jobfs_dir}"), "\n")

run_susie_gwas <- function(gene_name, ...) {
    loc <- df_gene_loc[ensembl_gene_id == gene_name]
    cis_start <- loc$cis_start
    cis_end <- loc$cis_end
    gwas_dt  <- df_gwas[position %between% c(cis_start, cis_end)]
    
    tmp_gwas    <- tempfile(tmpdir = jobfs_dir)
    tmp_ld_gwas <- tempfile(tmpdir = jobfs_dir)
    tmp_var_gwas <- tempfile(tmpdir = jobfs_dir)

    # Fallback cleanup on early return or error — catches anything missed below
    on.exit({
        unlink(tmp_gwas)
        unlink(tmp_var_gwas)
        unlink(list.files(dirname(tmp_ld_gwas),
                          pattern = paste0("^", basename(tmp_ld_gwas)),
                          full.names = TRUE))
    }, add = TRUE)

    writeLines(paste(gwas_dt$A1, gwas_dt$SNP, sep = "\t"), tmp_gwas)
    writeLines(gwas_dt$SNP, tmp_var_gwas)

    exit_code <- tryCatch(
        sys::exec_wait("plink", c(
            "--bfile",     paste0(dir_bfile, "/chr", chr),
            "--r",         "square",
            "--make-just-bim",
            "--chr",       chr,
            "--from-bp",   cis_start,
            "--to-bp",     cis_end,
            "--extract", tmp_var_gwas,
            "--silent",
            "--a1-allele", tmp_gwas, "1", "2",
            "--threads",   n_threads,
            "--out",       tmp_ld_gwas
        )),
        error = function(e) 1L
    )
    unlink(tmp_var_gwas)
    
    if (exit_code != 0 || !file.exists(paste0(tmp_ld_gwas, ".ld"))) return(NULL)

    bim_gwas    <- data.table::fread(paste0(tmp_ld_gwas, ".bim"), showProgress = FALSE)
    mat_ld_gwas <- as.matrix(data.table::fread(paste0(tmp_ld_gwas, ".ld"),
                                                col.names = bim_gwas$V2,
                                                showProgress = FALSE))
    rownames(mat_ld_gwas) <- bim_gwas$V2

    # Remove plink temp files immediately after mat_ld_gwas is in memory
    # on.exit() above is kept as a safety net for unexpected early returns
    unlink(tmp_gwas)
    unlink(list.files(dirname(tmp_ld_gwas),
                      pattern = paste0("^", basename(tmp_ld_gwas)),
                      full.names = TRUE))

    vars_na_gwas <- unique(rownames(which(is.na(mat_ld_gwas), arr.ind = TRUE)))
    vars <- intersect(bim_gwas$V2, gwas_dt$SNP)
    vars <- setdiff(vars, vars_na_gwas)
    mat_ld_gwas <- mat_ld_gwas[vars, vars]
    gwas_dt <- gwas_dt[SNP %in% vars][order(match(SNP, vars))]
    z <- gwas_dt[SNP %in% vars, b/se]
    N <- df_pheno_meta$n_eff[1]
    susie_res <- susie_rss(
        z = z,
        R = mat_ld_gwas,
        max_iter = opt$runsusie_maxit,
        n = N,
        # coloc set this as default
        ...
    )

    if (!susie_res$converged) return(NULL)
    susie_res <- annotate_susie(susie_res, vars, mat_ld_gwas)
    return(susie_res)
}

# run across all gwas genes and save results as RDS file
l_susie_gwas <- map(gwas_genes, function(g) {
    tryCatch(
        run_susie_gwas(g),
        error = function(e) {
            message(sprintf("[%s] ERROR gene=%s: %s", Sys.time(), g, e$message))
            NULL   # return NULL for this gene, continue to next
        }
    )
})
names(l_susie_gwas) <- gwas_genes

saveRDS(l_susie_gwas, file = output_file)
