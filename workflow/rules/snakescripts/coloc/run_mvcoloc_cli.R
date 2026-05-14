# CLI wrapper for nci-parallel
library(optparse)
library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(fst)
library(sys)
library(susieR)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--study",             type = "character"),
    make_option("--chr",               type = "character"),
    make_option("--pheno",             type = "character"),
    make_option("--dir_eqtl",          type = "character"),
    make_option("--ld_eqtl",           type = "character"),
    make_option("--susie_gwas",        type = "character"),
    make_option("--gene_loc",          type = "character"),
    make_option("--dir_bfile",         type = "character"),
    make_option("--output",            type = "character"),
    make_option("--threads",           type = "integer",  default = 8),
    make_option("--window_bp",         type = "integer",  default = 100000),
    make_option("--runsusie_coverage", type = "double",   default = 0.1),
    make_option("--p12",               type = "double",   default = 1e-5),
    make_option("--runsusie_maxit",    type = "integer",  default = 200),
    make_option("--runsusie_timeout",  type = "integer",  default = 180),
    make_option("--runsusie_repeat",   type = "logical",  default = FALSE),
    make_option("--coloc_timeout",     type = "integer",  default = 180)
)))

# Validate required args
required <- c("study", "chr", "pheno", "susie_gwas", "gene_loc", "dir_bfile", "output")
missing <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

cat(sprintf("[%s] PID %d | %s chr%s %s | threads=%d\n",
            Sys.time(), Sys.getpid(), opt$study, opt$chr, opt$pheno, opt$threads))

# Print all resolved paths for debugging
cat(sprintf("  susie_gwas:   %s (exists=%s)\n", opt$susie_gwas,    file.exists(opt$susie_gwas)))
cat(sprintf("  gene_loc:      %s (exists=%s)\n", opt$gene_loc,       file.exists(opt$gene_loc)))
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
susie_gwas_file   <- opt$susie_gwas
gene_loc_file     <- opt$gene_loc
dir_bfile         <- opt$dir_bfile
output_file       <- opt$output
runsusie_coverage <- opt$runsusie_coverage
runsusie_maxit    <- opt$runsusie_maxit
runsusie_repeat   <- opt$runsusie_repeat
coloc_p12 <- opt$p12

dir_eqtl <- opt$dir_eqtl
ld_eqtl <- opt$ld_eqtl

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
    formatted <- paste0("[", Sys.time(), "] [PID:", Sys.getpid(), "] ", msg)
    cat(formatted, "\n")
    tryCatch(
        write(formatted, file = progress_log, append = TRUE),
        error = function(e) invisible(NULL)
    )
}

log_progress(glue("Start: study={study} chr={chr} pheno={pheno} threads={n_threads}"))

io_threads <- min(n_threads, 8)
setDTthreads(io_threads)

cat(glue("Loading data with {io_threads} threads..."), "\n")

# --- load SuSiE GWAS files (created in previous step) ---
t0 <- Sys.time()

l_susie_gwas_full <- readRDS(susie_gwas_file)

# load eQTL files
eqtl_files <- list.files(
    path    = opt$dir_eqtl,
    pattern = glue("chr{opt$chr}\\.fst$"),
    full.names  = TRUE,
    recursive   = TRUE
)

if (length(eqtl_files) == 0) {
    log_progress(glue("No eQTL .fst files found for chr{opt$chr} in {opt$dir_eqtl} — writing empty output"))
    fwrite(data.table(), opt$output, sep = "\t")
    quit(save = "no", status = 0)
}
names(eqtl_files) <- fs::path_split(eqtl_files) |> map_chr(4)

df_eqtl <- map_df(eqtl_files, read_fst, as.data.table = TRUE, .id = "biosample") |>
    filter(gene %in% names(l_susie_gwas_full))
setkey(df_eqtl, gene, biosample)

# Keep only genes that are in both GWAS susie AND eQTL data
egene_list   <- intersect(names(l_susie_gwas_full), unique(df_eqtl$gene))

l_susie_gwas <- l_susie_gwas_full[egene_list]   # subset to needed genes only
rm(l_susie_gwas_full)                            # free immediately
gc(verbose = FALSE)

cat(glue("  Retained {length(egene_list)} genes present in both GWAS SuSiE and eQTL data\n"))


# --- Gene locations ---
t0 <- Sys.time()
CHR <- chr
df_gene_loc <- fread(gene_loc_file)[chr == CHR & ensembl_gene_id %in% egene_list]
df_gene_loc[, `:=`(cis_start = pmax(1L, start - window_bp), cis_end = end + window_bp)]
cat(glue("  Gene locations: {nrow(df_gene_loc)} genes ({round(difftime(Sys.time(), t0, units='secs'), 1)}s)"), "\n")

# ==========================================================================
# 3. Run Coloc SuSiE
# ==========================================================================
jobfs_dir <- Sys.getenv("PBS_JOBFS", tempdir())
cat(glue("  jobfs_dir: {jobfs_dir}"), "\n")

safe_runsusie <- function(dc) {
    tryCatch(
        coloc::runsusie(dc,
            coverage                 = runsusie_coverage,
            maxit                    = runsusie_maxit,
            repeat_until_convergence = runsusie_repeat,
            estimate_prior_variance  = FALSE),
        error = function(e) NULL
    )
}

safe_runsusie <- function(dc, timeout_secs = opt$runsusie_timeout) {
    # setTimeLimit raises an error if the function exceeds timeout_secs
    # This catches infinite loops in susieR that tryCatch(error=) alone cannot catch
    setTimeLimit(cpu = timeout_secs, elapsed = timeout_secs, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)

    tryCatch(
        coloc::runsusie(dc,
            coverage                 = runsusie_coverage,
            maxit                    = runsusie_maxit,
            repeat_until_convergence = runsusie_repeat,
            estimate_prior_variance  = FALSE),
        error = function(e) {
            if (grepl("time limit|reached elapsed", conditionMessage(e), ignore.case = TRUE)) {
                log_progress(glue("  safe_runsusie TIMEOUT after {timeout_secs}s — skipping"))
            } else {
                log_progress(glue("  safe_runsusie ERROR: {conditionMessage(e)}"))
            }
            NULL
        }
    )
}

safe_coloc_susie <- function(susie_eqtl, susie_gwas, p12, timeout_secs = opt$coloc_timeout) {
    # coloc.susie can also hang if susie objects have degenerate credible sets
    setTimeLimit(cpu = timeout_secs, elapsed = timeout_secs, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)

    tryCatch(
        coloc::coloc.susie(
            dataset1 = susie_eqtl,
            dataset2 = susie_gwas,
            p12      = p12
        ),
        error = function(e) {
            if (grepl("time limit|reached elapsed", conditionMessage(e), ignore.case = TRUE)) {
                log_progress(glue("  coloc.susie TIMEOUT after {timeout_secs}s — skipping"))
            } else {
                log_progress(glue("  coloc.susie ERROR: {conditionMessage(e)}"))
            }
            NULL
        }
    )
}

load_ld_gene <- function(gene) {
    bim_file <- paste0(ld_eqtl, "/", gene, ".bim")
    ld_file  <- paste0(ld_eqtl, "/", gene, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) return(NULL)

    bim_eqtl    <- data.table::fread(bim_file, showProgress = FALSE)
    mat_ld_eqtl <- as.matrix(data.table::fread(ld_file,
                                                col.names = bim_eqtl$V2,
                                                showProgress = FALSE))
    rownames(mat_ld_eqtl) <- bim_eqtl$V2

    # Remove missing variants 
    vars_na <- unique(rownames(which(is.na(mat_ld_eqtl), arr.ind = TRUE)))
    if (length(vars_na) > 0) {
        keep <- setdiff(rownames(mat_ld_eqtl), vars_na)
        mat_ld_eqtl <- mat_ld_eqtl[keep, keep, drop = FALSE]
    }

    if (nrow(mat_ld_eqtl) < 10) return(NULL)

    return(mat_ld_eqtl)
}

df_to_dc <- function(df, vars){
    df <- filter(df, snp %in% vars) |> 
        arrange(match(snp, vars))
    list(
        beta     = df$beta,
        varbeta  = df$varbeta,
        snp      = df$snp,
        position = df$position,
        MAF      = df$MAF,
        type     = "quant",
        N        = df$N[1]
    )
}

# REMOVE: genes_ld preloaded for all genes
# genes_ld <- map(egene_list, load_ld_gene)
# names(genes_ld) <- egene_list

log_progress(glue("Processing {length(egene_list)} genes sequentially"))
setDTthreads(n_threads)

all_results <- vector("list", length(egene_list))
names(all_results) <- egene_list

for (i in seq_along(egene_list)) {
    gene_name <- egene_list[i]

    # Declare cleanup at the top — runs on every exit path from this iteration
    # Use a local() block to scope on.exit() to this iteration only
    local({
        mat_ld_eqtl <- NULL
        susie_gwas  <- NULL
        eqtl_gene   <- NULL

        on.exit({
            rm(mat_ld_eqtl, susie_gwas, eqtl_gene)
            gc(verbose = FALSE)
        }, add = TRUE)

        mat_ld_eqtl <- load_ld_gene(gene_name)
        if (is.null(mat_ld_eqtl)) return(invisible(NULL))   # next → return in local()

        susie_gwas <- l_susie_gwas[[gene_name]]
        if (is.null(susie_gwas)) return(invisible(NULL))

        l_susie_gwas[[gene_name]] <<- NULL   # <<- to modify outer environment

        eqtl_gene <- df_eqtl[gene == gene_name]
        bs_names  <- unique(eqtl_gene$biosample)
        bs_results <- vector("list", length(bs_names))

        for (j in seq_along(bs_names)) {
            bs      <- bs_names[j]
            eqtl_dt <- eqtl_gene[biosample == bs]

            vars <- intersect(eqtl_dt$snp, rownames(mat_ld_eqtl))
            if (length(vars) < 10) next

            dc     <- df_to_dc(eqtl_dt, vars)
            dc$LD  <- mat_ld_eqtl[vars, vars, drop = FALSE]

            t_bs       <- Sys.time()
            susie_eqtl <- safe_runsusie(dc)
            rm(dc)

            if (is.null(susie_eqtl)) {
                log_progress(glue("  [{i}/{length(egene_list)}] gene={gene_name} bs={bs} | susie_eqtl=NULL, skipping"))
                next
            }

            res     <- safe_coloc_susie(susie_eqtl, susie_gwas, coloc_p12)
            elapsed <- round(difftime(Sys.time(), t_bs, units = "secs"), 1)
            rm(susie_eqtl)

            if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0) {
                log_progress(glue("  [{i}/{length(egene_list)}] gene={gene_name} bs={bs} | no results ({elapsed}s)"))
                next
            }

            n_hits <- sum(res$summary$PP.H4.abf > 0.8, na.rm = TRUE)
            log_progress(glue("  [{i}/{length(egene_list)}] gene={gene_name} bs={bs} | {nrow(res$summary)} signals, PP.H4>0.8: {n_hits} ({elapsed}s)"))

            bs_results[[j]] <- cbind(
                data.table(biosample = bs, gene = gene_name),
                as.data.table(res$summary)
            )
        }

        bs_results <- Filter(Negate(is.null), bs_results)
        if (length(bs_results) > 0) {
            all_results[[i]] <<- rbindlist(bs_results, use.names = TRUE, fill = TRUE)
        }
    })  # on.exit() fires here — all paths guaranteed

    if (i %% 50 == 0 || i == length(egene_list)) {
        n_done <- sum(!sapply(all_results, is.null))
        log_progress(glue("Progress: {i}/{length(egene_list)} genes, {n_done} with results"))
    }
}

df_results <- rbindlist(Filter(Negate(is.null), all_results),
                        use.names = TRUE, fill = TRUE)

fwrite(df_results, output_file, sep = "\t")
log_progress(glue("Done: {nrow(df_results)} rows written to {output_file}"))
