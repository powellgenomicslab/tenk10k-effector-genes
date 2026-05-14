# Multivariant coloc analysis with SuSiE
# Supports both Snakemake S4 object (via script: directive) and
# plain list (via CLI wrapper for nci-parallel)

library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(fs)
library(fst)
library(sys)
library(susieR)

# ==========================================================================
# Normalise snakemake object: S4 (@) → plain list ($)
# ==========================================================================
if (isClass("Snakemake") && is(snakemake, "Snakemake")) {
    smk <- list(
        input     = setNames(as.list(snakemake@input), names(snakemake@input)),
        output    = setNames(as.list(snakemake@output), names(snakemake@output)),
        params    = setNames(as.list(snakemake@params), names(snakemake@params)),
        wildcards = setNames(as.list(snakemake@wildcards), names(snakemake@wildcards)),
        threads   = snakemake@threads
    )
} else {
    smk <- snakemake
}

# ==========================================================================
# 1. Parameters
# ==========================================================================
dir_eqtl          <- smk$input[["dir_eqtl"]]
ld_eqtl           <- smk$input[["ld_eqtl"]]
gwas_file         <- smk$input[["gwas"]]
gene_loc_file     <- smk$input[["gene_loc"]]
pheno_metadata_file <- smk$input[["pheno_metadata"]]
dir_bfile         <- smk$input[["dir_bfile"]]
output_file       <- smk$output[["coloc"]]

PARAMS   <- smk$params
window_bp <- PARAMS[["window_bp"]]
n_threads <- smk$threads

pheno <- smk$wildcards[["pheno"]]
study <- smk$wildcards[["study"]]
chr   <- smk$wildcards[["chr"]]

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

# ==========================================================================
# 2. Load lightweight data only — NO eQTL rows loaded yet
# ==========================================================================
io_threads <- min(n_threads, 8)
setDTthreads(io_threads)
fst::threads_fst(nr_of_threads = min(io_threads, 8), reset_after_fork = FALSE)

cat(glue("Loading data with {io_threads} threads..."), "\n")

# --- Phenotype metadata ---
t0 <- Sys.time()
df_pheno_meta <- fread(pheno_metadata_file)[include == TRUE & trait_id == pheno]
pheno_type <- ifelse(df_pheno_meta$supercategory == "biological", "quant", "cc")
pheno_sd   <- ifelse(pheno_type == "quant", 1, NA)
cat(glue("  Pheno metadata: {round(difftime(Sys.time(), t0, units='secs'), 1)}s"), "\n")

# --- eQTL file index — paths only, NO data loaded ---
t0 <- Sys.time()
eqtl_files <- dir_ls(dir_eqtl, recurse = TRUE, glob = glue("*chr{chr}.fst"))
eqtl_files <- set_names(
    eqtl_files,
    eqtl_files |> str_remove_all(dir_eqtl) |> path_split() |> map_chr(2)
)

df_eqtl <- map_df(eqtl_files, ~read_fst(.x, as.data.table = TRUE, ), .id = "biosample")
setkey(df_eqtl, snp)

egene_list <- unique(df_eqtl$gene)
cat(glue("  Index: {length(egene_list)} genes across {length(eqtl_files)} biosamples ",
         "({round(difftime(Sys.time(), t0, units='secs'), 1)}s)"), "\n")

# --- Gene locations ---
t0 <- Sys.time()
CHR <- chr
df_gene_loc <- fread(gene_loc_file)[chr == CHR]
df_gene_loc[, `:=`(cis_start = pmax(1L, start - window_bp), cis_end = end + window_bp)]
setkey(df_gene_loc, ensembl_gene_id)
cat(glue("  Gene locations: {nrow(df_gene_loc)} genes ({round(difftime(Sys.time(), t0, units='secs'), 1)}s)"), "\n")

# --- GWAS — load full chr, keep SNPs seen in any eQTL file ---
t0 <- Sys.time()
df_gwas <- fread(gwas_file)[SNP %in% unique(df_eqtl$snp)]
dup_snps <- df_gwas[duplicated(SNP), unique(SNP)]
df_gwas <- df_gwas[!SNP %in% dup_snps]
df_gwas[, `:=`(varbeta = se^2, MAF = pmin(freq, 1 - freq))]
setnames(df_gwas, "SNP", "snp")
setnames(df_gwas, "b",   "beta")
setkey(df_gwas, snp)
df_gwas[df_eqtl, position := i.position, on = "snp"]

cat(glue("  GWAS: {nrow(df_gwas)} rows ({round(difftime(Sys.time(), t0, units='secs'), 1)}s)"), "\n")

# filter only region with egene and P_gwas < threshold
df_gwas_top_pos <- df_gwas[P < PARAMS$min_p_gwas, .(position)]
matched_indices <- df_gene_loc[df_gwas_top_pos, on = .(cis_start <= position, cis_end >= position), which = TRUE, nomatch = 0]
gwas_genes <- unique(df_gene_loc[matched_indices])$ensembl_gene_id

df_eqtl <- df_eqtl[df_gwas[,.(snp)], nomatch = 0L][gene %in% gwas_genes]

egene_list <- unique(df_eqtl$gene)

rm(dup_snps)
gc()

log_progress(glue("Loaded: {length(egene_list)}, {nrow(df_gwas)} GWAS SNPs, {nrow(df_eqtl)} eQTL SNPS * 28 cell types"))
cat(glue("  Base RSS: eQTL = {format(object.size(df_eqtl), units='auto')} +",
         " GWAS ={format(object.size(df_gwas), units='auto')}\n\n"))

# ==========================================================================
# 3. No pre-split — eQTL loaded on demand inside run_mv_coloc()
# ==========================================================================
jobfs_dir <- Sys.getenv("PBS_JOBFS", tempdir())
cat(glue("  jobfs_dir: {jobfs_dir}"), "\n")

# ==========================================================================
# 6. Per-gene coloc — loads eQTL from fst ON DEMAND, frees after each gene
# ==========================================================================
safe_runsusie <- function(z, LD) {
    tryCatch(
        coloc::runsusie(dc,
            coverage                 = runsusie_coverage,
            maxit                    = runsusie_maxit,
            repeat_until_convergence = runsusie_repeat,
            estimate_prior_variance  = FALSE),
        error = function(e) NULL
    )
}

# run susie on GWAS datasets for all possible eGenes

l_susie_gwas <- map(egene_list, ~ {
    gene_name <- .x
    loc <- df_gene_loc[ensembl_gene_id == gene_name]
    cis_start <- loc$cis_start
    cis_end <- loc$cis_end
    gwas_dt  <- df_gwas[position %between% c(cis_start, cis_end)]
    
    # ------------------------------------------------------------------
    # Compute GWAS LD on demand via plink
    # ------------------------------------------------------------------
    tmp_gwas    <- tempfile(tmpdir = jobfs_dir)
    tmp_ld_gwas <- tempfile(tmpdir = jobfs_dir)
    on.exit({
        unlink(tmp_gwas)
        unlink(list.files(dirname(tmp_ld_gwas),
                          pattern = basename(tmp_ld_gwas),
                          full.names = TRUE))
    }, add = TRUE)

    writeLines(paste(gwas_dt$A1, gwas_dt$snp, sep = "\t"), tmp_gwas)

    exit_code <- tryCatch(
        sys::exec_wait("plink", c(
            "--bfile",     paste0(dir_bfile, "/chr", chr),
            "--r",         "square",
            "--make-just-bim",
            "--chr",       chr,
            "--from-bp",   cis_start,
            "--to-bp",     cis_end,
            "--silent",
            "--a1-allele", tmp_gwas, "1", "2",
            "--threads",   n_threads,
            "--out",       tmp_ld_gwas
        )),
        error = function(e) 1L
    )

    if (exit_code != 0 || !file.exists(paste0(tmp_ld_gwas, ".ld"))) return(NULL)

    bim_gwas    <- data.table::fread(paste0(tmp_ld_gwas, ".bim"), showProgress = FALSE)
    mat_ld_gwas <- as.matrix(data.table::fread(paste0(tmp_ld_gwas, ".ld"),
                                                col.names = bim_gwas$V2,
                                                showProgress = FALSE))
    rownames(mat_ld_gwas) <- bim_gwas$V2

    vars_na_gwas <- unique(rownames(which(is.na(mat_ld_gwas), arr.ind = TRUE)))
    vars <- intersect(bim_gwas$V2, gwas_dt$snp)
    vars <- setdiff(vars, vars_na_gwas)
    mat_ld_gwas <- mat_ld_gwas[vars, vars]
    gwas_dt <- gwas_dt[snp %in% vars][order(match(snp, vars))]
    z <- gwas_dt[snp %in% vars, beta/se]
    N <- df_pheno_meta$n_eff[1]
    susie_res <- susie_rss(z = z, R = mat_ld_gwas, max_iter = PARAMS$runsusie_maxit, n = N)

    if (!susie_res$converged) return(NULL)
    susie_res <- annotate_susie(susie_res, vars, mat_ld_gwas)
})

run_mv_coloc <- function(gene_name, dt_threads = n_threads, plink_threads = n_threads) {
    setDTthreads(dt_threads)

    eqtl_dt <- df_eqtl[gene == gene_name]

    # ------------------------------------------------------------------
    # Subset GWAS to cis-window (reference df_gwas — no copy of full table)
    # ------------------------------------------------------------------
    cis_start <- eqtl_dt[,min(position)]
    cis_end   <- eqtl_dt[,max(position)]
    gwas_dt  <- df_gwas[position %between% c(cis_start, cis_end)]
    if (nrow(gwas_dt) == 0) return(NULL)

    # ------------------------------------------------------------------
    # Load eQTL LD on demand
    # ------------------------------------------------------------------
    bim_file <- paste0(ld_eqtl, "/", gene_name, ".bim")
    ld_file  <- paste0(ld_eqtl, "/", gene_name, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) return(NULL)

    bim_eqtl    <- data.table::fread(bim_file, showProgress = FALSE)
    mat_ld_eqtl <- as.matrix(data.table::fread(ld_file,
                                                col.names = bim_eqtl$V2,
                                                showProgress = FALSE))
    rownames(mat_ld_eqtl) <- bim_eqtl$V2

    # ------------------------------------------------------------------
    # Compute GWAS LD on demand via plink
    # ------------------------------------------------------------------
    tmp_gwas    <- tempfile(tmpdir = jobfs_dir)
    tmp_ld_gwas <- tempfile(tmpdir = jobfs_dir)
    on.exit({
        unlink(tmp_gwas)
        unlink(list.files(dirname(tmp_ld_gwas),
                          pattern = basename(tmp_ld_gwas),
                          full.names = TRUE))
    }, add = TRUE)

    writeLines(paste(gwas_dt$A1, gwas_dt$snp, sep = "\t"), tmp_gwas)

    exit_code <- tryCatch(
        sys::exec_wait("plink", c(
            "--bfile",     paste0(dir_bfile, "/chr", chr),
            "--r",         "square",
            "--make-just-bim",
            "--chr",       chr,
            "--from-bp",   cis_start,
            "--to-bp",     cis_end,
            "--silent",
            "--a1-allele", tmp_gwas, "1", "2",
            "--threads",   plink_threads,
            "--out",       tmp_ld_gwas
        )),
        error = function(e) 1L
    )

    if (exit_code != 0 || !file.exists(paste0(tmp_ld_gwas, ".ld"))) return(NULL)

    bim_gwas    <- data.table::fread(paste0(tmp_ld_gwas, ".bim"), showProgress = FALSE)
    mat_ld_gwas <- as.matrix(data.table::fread(paste0(tmp_ld_gwas, ".ld"),
                                                col.names = bim_gwas$V2,
                                                showProgress = FALSE))
    rownames(mat_ld_gwas) <- bim_gwas$V2

    # ------------------------------------------------------------------
    # Variant overlap
    # ------------------------------------------------------------------
    vars_na_eqtl <- unique(rownames(which(is.na(mat_ld_eqtl), arr.ind = TRUE)))
    vars_na_gwas <- unique(rownames(which(is.na(mat_ld_gwas), arr.ind = TRUE)))
    vars <- Reduce(intersect, list(eqtl_dt$snp, gwas_dt$snp, bim_eqtl$V2, bim_gwas$V2))
    vars <- setdiff(vars, c(vars_na_eqtl, vars_na_gwas))
    if (length(vars) < 10) return(NULL)

    # ------------------------------------------------------------------
    # GWAS SuSiE (computed once, shared across biosamples)
    # ------------------------------------------------------------------
    gwas_dt  <- gwas_dt[snp %in% vars]
    
    data_gwas <- gwas_dt[,.(beta, varbeta, snp, position, MAF)] |> 
     as.list()

    data_gwas$type <- pheno_type
    data_gwas$LD <- mat_ld_gwas[vars, vars]
    data_gwas$N <- df_pheno_meta$n_eff
    if (pheno_type == "quant") data_gwas$sdY <- pheno_sd

    rm(mat_ld_gwas, bim_gwas, gwas_dt)

    susie_gwas <- safe_runsusie(data_gwas)
    
    if (is.null(susie_gwas)) return(NULL)
    rm(data_gwas)

    biosamples <- unique(eqtl$biosample)

    l_eqtl <- df_eqtl[gene == gene_name, .(data = list(.SD), by = biosample)] |> 
        mutate(susie = map(data, ~ {
                    dc <- as.list(.x[, .(beta, varbeta, snp, position, MAF)])
                    dc$type <- "quant"
                    dc$ld <- mat_ld_eqtl[dc$snp, dc$snp, drop = FALSE]
                    dc$N <- .x$N[1]
                    safe_runsusie(dc)
                }))
    results <- map(biosamples, function(bs) {
        eqtl  <- which(df_eqtl$biosample == bs & eqtl$snp %in% vars)
        snps <- eqtl$snp[idx]

        data_eqtl <- df_eqtl[biosample == bs & gene == gene_name] |> 

            list(
            beta     = eqtl$beta[idx],
            varbeta  = eqtl$varbeta[idx],
            snp      = snps,
            position = eqtl$position[idx],
            MAF      = eqtl$MAF[idx],
            type     = "quant",
            LD       = mat_ld_eqtl[snps, snps, drop = FALSE],
            N        = eqtl$N[idx[1]]
        )

        susie_eqtl <- safe_runsusie(data_eqtl)
        rm(data_eqtl)
        if (is.null(susie_eqtl)) return(NULL)

        t_coloc <- Sys.time()
        cat(glue("    [{format(t_coloc, '%H:%M:%S')}] coloc.susie: gene={gene_name} bs={bs}"), "\n")
        res <- tryCatch(
            coloc::coloc.susie(
            dataset1 = susie_eqtl,
            dataset2 = susie_gwas,
            p12      = p12_param
            ),
            error = function(e) NULL
        )
        rm(susie_eqtl)

        if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0) {
            cat(glue("    [{format(Sys.time(), '%H:%M:%S')}] coloc.susie null/empty result: gene={gene_name} bs={bs}"), "\n")
            return(NULL)
        }
        cbind(
            data.table::data.table(
                chr = as.integer(chr), biosample = bs,
                pheno = pheno, gene = gene_name),
            data.table::as.data.table(res$summary)
        )

        cat(glue("    [{format(Sys.time(), '%H:%M:%S')}] coloc.susie done: {round(difftime(Sys.time(), t_coloc, units='secs'), 1)}s"), "\n")

    })

    rm(mat_ld_eqtl, bim_eqtl, susie_gwas, eqtl)

    results <- Filter(Negate(is.null), results)
    if (length(results) == 0) return(NULL)
    data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
}

# ==========================================================================
# 7. Process all genes SEQUENTIALLY — nci-parallel handles parallelism
# ==========================================================================
runsusie_coverage <- PARAMS[["runsusie_coverage"]]
runsusie_maxit    <- PARAMS[["runsusie_maxit"]]
runsusie_repeat   <- PARAMS[["runsusie_repeat"]]
p12_param         <- PARAMS[["p12"]]

log_progress(paste0("Processing ", length(egene_list), " genes sequentially (1 thread)"))

all_results <- list()
gene_status <- setNames(rep("pending", length(egene_list)), egene_list)

for (i in seq_along(egene_list)) {
    g <- egene_list[i]

    if (i %% 50 == 0 || i == length(egene_list)) {
        n_done  <- sum(gene_status %in% c("done", "done_empty"))
        n_error <- sum(gene_status == "error")
        cat(glue("  Gene {i}/{length(egene_list)} ({g}): {n_done} done, {n_error} errors"), "\n")
        flush(stdout())
    }

    res <- tryCatch(
        run_mv_coloc(g),
        error = function(e) {
            log_progress(glue("ERROR gene={g}: {conditionMessage(e)}"))
            structure(conditionMessage(e), class = "worker-error")
        }
    )

    if (inherits(res, "worker-error")) {
        gene_status[g] <- "error"
    } else if (is.null(res)) {
        gene_status[g] <- "done_empty"
    } else {
        gene_status[g] <- "done"
        all_results[[g]] <- res
    }

    # Periodic gc to prevent memory creep
    if (i %% 100 == 0) gc(verbose = FALSE)
}

n_done    <- sum(gene_status %in% c("done", "done_empty"))
n_error   <- sum(gene_status == "error")
n_with_res <- sum(gene_status == "done")
log_progress(glue("Completed: {n_done}/{length(egene_list)} genes, {n_error} errors, {n_with_res} with results"))

# ==========================================================================
# 8. Gene-level summary and write results
# ==========================================================================
gene_summary <- data.table::data.table(
    gene   = egene_list,
    status = gene_status[egene_list]
)

summary_tab <- table(gene_summary$status)
cat(glue("\n  Gene summary: {paste(names(summary_tab), summary_tab, sep='=', collapse=', ')}"), "\n")
log_progress(paste0("Gene summary: ", paste(names(summary_tab), summary_tab, sep = "=", collapse = ", ")))

# Warn if unexpected errors occurred
n_errors <- sum(gene_status == "error")
if (n_errors > 0) {
    warning(glue("{n_errors} genes failed with errors — check progress log: {progress_log}"))
}

# Write gene status log
status_file <- paste0(output_file, ".gene_status.tsv")
data.table::fwrite(gene_summary, status_file, sep = "\t")
cat(glue("  Gene status: {status_file}"), "\n")

# Combine results — only genes with actual coloc hits
setDTthreads(n_threads)
dir_create(dirname(output_file))

if (length(all_results) > 0) {
    df_coloc_res <- data.table::rbindlist(all_results, use.names = TRUE, fill = TRUE)
    cat(glue("  Results: {nrow(df_coloc_res)} rows from {length(all_results)} genes"), "\n")
} else {
    # Header-only file so downstream knows the task completed
    df_coloc_res <- data.table::data.table(
        chr = integer(0), biosample = character(0), pheno = character(0),
        gene = character(0), nsnps = integer(0), hit1 = character(0),
        hit2 = character(0), PP.H0.abf = numeric(0), PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0), PP.H3.abf = numeric(0), PP.H4.abf = numeric(0),
        idx1 = integer(0), idx2 = integer(0)
    )
    cat("  Results: 0 rows (no coloc hits for any gene)\n")
}

data.table::fwrite(df_coloc_res, output_file, row.names = FALSE, sep = "\t")
cat(glue("  Written: {output_file} ({nrow(df_coloc_res)} rows)\n"))
log_progress(glue("Results written: {output_file} ({nrow(df_coloc_res)} rows)"))