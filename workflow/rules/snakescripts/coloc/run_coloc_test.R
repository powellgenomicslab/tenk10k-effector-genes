# Coloc analysis script for Snakemake pipeline
# Modified from original run_coloc.R to work with snakemake

# Load libraries
library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(fs)
library(fst)

## 1. Get Snakemake parameters -------------------------------------------------------------
# Input files
eqtl_file <- snakemake@input[["eqtl"]]
gwas_file <- snakemake@input[["gwas"]]
gene_loc_file <- snakemake@input[["gene_loc"]]
pheno_metadata_file <- snakemake@input[["pheno_metadata"]]

# Output file
output_file <- snakemake@output[["coloc"]]

# Logging. Without this nothing reaches the `log:` path declared in coloc.smk,
# which is why logs/coloc/ did not exist at all and per-gene coloc errors left no
# trace. Adapted from snakescripts/enrichment/{stringdb,gprofiler}.R but WITHOUT
# their on.exit() call - see the note below; those two scripts have the same
# latent bug and their logs are likely empty for the same reason.
log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, append = TRUE, type = "message")
# NB: do NOT use on.exit() here. snakemake source()s this script, and source()
# eval()s one expression at a time, so on.exit() attaches to that short-lived
# eval frame and fires immediately - closing the log before anything is written
# (the reason logs/coloc/ was empty). Call close_log() explicitly instead, on
# every exit path. On an uncaught error the sink stays open and the error text
# lands in the log, which is what we want; R flushes connections on exit.
.log_closed <- FALSE
close_log <- function() {
    if (.log_closed) return(invisible(NULL))
    .log_closed <<- TRUE
    sink(type = "message")
    sink(type = "output")
    close(log_con)
}

# Parameters
window_bp <- snakemake@params[["window_bp"]]

# Extract metadata from wildcards
biosample <- snakemake@wildcards[["biosample"]]
pheno <- snakemake@wildcards[["pheno"]]
study <- snakemake@wildcards[["study"]]
chr <- snakemake@wildcards[["chr"]]

# phenotype metadata.
#
# trait_metadata_n.tsv intentionally holds more than one row per trait_id where a
# trait has multiple source studies (the BMD traits carry both Kemp2024 and
# Miao2024); the intended one is flagged `include == TRUE`. This script was
# missing that filter, unlike its siblings run_susie_gwas_cli.R and
# run_multivariant_coloc.R, so `ifelse(df_pheno_meta$supercategory == ...)`
# returned a length-2 vector and `if (pheno_type == "quant")` below aborted the
# whole script with "the condition has length > 1" (an error since R 4.2, only a
# warning before, which is why these traits ran under older R).
df_pheno_meta <- fread(pheno_metadata_file) %>%
    filter(include == TRUE, trait_id == pheno)
if (nrow(df_pheno_meta) != 1) {
    stop(sprintf("expected exactly 1 included trait_metadata row for %s, got %d",
                 pheno, nrow(df_pheno_meta)))
}
pheno_type <- if (df_pheno_meta$supercategory == "biological") "quant" else "cc"

# assume that quantitative traits have sdY = 1
# todo: add sdY in the metadata file if needed
pheno_sd <- ifelse(pheno_type == "quant", 1, NA)

## 2. Load data -------------------------------------------------------------
df_eqtl <- read_fst(eqtl_file, as.data.table = TRUE)
setkey(df_eqtl, snp)
egene_list <- unique(df_eqtl$gene)

df_gwas <- fread(gwas_file)
setkey(df_gwas, SNP)

# Clean the GWAS before joining. Two data-quality issues in the .ma files each
# made coloc.abf fail for *every* gene on a chromosome, so the chromosome
# vanished from the output entirely:
#
#  1. b = 0 with se = 0, used by some GWAS for unreported effects (sle: 48,700
#     of 581,834 SNPs on chr1). z = 0/0 = NaN, so min(p) is NA and coloc's
#     check_dataset() dies on `if (min(p) > warn.minp)` with the opaque
#     "missing value where TRUE/FALSE needed".
#  2. Duplicate SNP IDs (pancreatic: 25,833 on chr3). The keyed join below
#     emits one row per duplicate, multiplying the eQTL rows, and
#     check_dataset() then rejects dataset 1 with "duplicated snps found".
#
# One bad SNP anywhere in a cis window was enough to discard that gene, which
# is why affected traits lost whole chromosomes while clean traits (asthma) lost
# none.
#
# Duplicate IDs are dropped outright rather than de-duplicated to the first
# occurrence: most carry *conflicting* records, i.e. distinct variants collapsed
# onto one ID, each with its own freq/b/se (ra chr3: 6,720 of 6,723 duplicate
# IDs conflict; parkinsons: 36,729 of 36,831). Keeping an arbitrary one would
# silently attribute one variant's effect to the shared ID, so an ambiguous ID
# contributes nothing.
n_gwas_all <- nrow(df_gwas)
df_gwas <- df_gwas[is.finite(b) & is.finite(se) & se > 0]
n_bad_se <- n_gwas_all - nrow(df_gwas)
df_gwas <- df_gwas[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
n_dup <- n_gwas_all - n_bad_se - nrow(df_gwas)
if (n_bad_se > 0 || n_dup > 0) {
    cat(glue("GWAS cleaning: dropped {n_bad_se} SNPs with non-finite/zero se and ",
             "{n_dup} rows with ambiguous (duplicated) SNP IDs, from {n_gwas_all}"), "\n")
}
setkey(df_gwas, SNP)

# eQTL side: guard against duplicate gene-SNP pairs for the same reason.
# unique() with a `by` other than the key drops the key, so re-set it - the
# joins below rely on df_eqtl being keyed on snp.
df_eqtl <- unique(df_eqtl, by = c("gene", "snp"))
setkey(df_eqtl, snp)

# subset gwas & eqtl data
df_gwas_subset <- df_gwas[df_eqtl, nomatch=NULL] %>%
    mutate(varbeta = se^2) %>%
    select(gene, beta = b, varbeta, snp = SNP, position)

df_eqtl_subset <- df_eqtl[df_gwas[,.(SNP)], nomatch=NULL]

df_gene_loc <- fread(gene_loc_file) %>%
    filter(ensembl_gene_id %in% egene_list) %>%
    mutate(
        cis_start = pmax(1, start - window_bp),
        cis_end = end + window_bp
    )

# Single output schema for every code path. `status` is "ok" for a successful
# test, otherwise the coloc.abf error message.
#
# A failed gene emits a row (df_failed) rather than vanishing. Previously the
# error path returned df_blank(), which has zero rows, so a chromosome where
# every gene errored was indistinguishable from one with no eligible eGenes -
# both wrote a header-only file and the chromosome silently disappeared from
# all_chr.coloc.tsv.
df_blank <- function(){
    data.table(
        chr = character(0),
        biosample = character(0),
        pheno = character(0),
        gene = character(0),
        nsnps_coloc_tested = numeric(0),
        PP.H0.abf = numeric(0),
        PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0),
        PP.H3.abf = numeric(0),
        PP.H4.abf = numeric(0),
        top_snp = character(0),
        top_snp_pph4 = numeric(0),
        status = character(0)
    )
}

df_failed <- function(gene_name, msg){
    data.table(
        chr = chr,
        biosample = biosample,
        pheno = pheno,
        gene = gene_name,
        nsnps_coloc_tested = NA_real_,
        PP.H0.abf = NA_real_,
        PP.H1.abf = NA_real_,
        PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_,
        PP.H4.abf = NA_real_,
        top_snp = NA_character_,
        top_snp_pph4 = NA_real_,
        status = msg
    )
}

if (nrow(df_gene_loc) == 0) {
    cat(glue("No eGenes found for chromosome {chr} - creating empty output"), "\n")
    # Create empty output with correct structure
    dir_create(dirname(output_file))
    fwrite(df_blank(), output_file, row.names = FALSE, sep = "\t")
    close_log()   # flush the sink before quitting on this early-exit path
    quit(save = "no")
}

## 3. Colocalisation analysis -------------------------------------------------------------
run_coloc <- function(gene_name) {
    cis_start <- df_gene_loc[ensembl_gene_id == gene_name, cis_start]
    cis_end <- df_gene_loc[ensembl_gene_id == gene_name, cis_end]
    df_eqtl_gene <- df_eqtl_subset[gene == gene_name & position %between% c(cis_start, cis_end)]
    df_gwas_gene <- df_gwas_subset[gene == gene_name & position %between% c(cis_start, cis_end)]
    data_eqtl <- as.list(df_eqtl_gene[, !"gene"])
    data_eqtl$type <- "quant"
    data_gwas <- as.list(df_gwas_gene[, !"gene"])
    data_gwas$type <- pheno_type

    if (pheno_type == "quant") data_gwas$sdY <- pheno_sd

    # Perform colocalisation analysis
    tryCatch({
        my.res <- coloc.abf(
            dataset1 = data_eqtl,
            dataset2 = data_gwas
        )
        
        # Extract results
       data.table(
            chr = chr,
            biosample = biosample,
            pheno = pheno,
            gene = gene_name,
            nsnps_coloc_tested = my.res$summary[1],
            PP.H0.abf = my.res$summary[2],
            PP.H1.abf = my.res$summary[3],
            PP.H2.abf = my.res$summary[4],
            PP.H3.abf = my.res$summary[5],
            PP.H4.abf = my.res$summary[6],
            top_snp = my.res$results[which.max(my.res$results$SNP.PP.H4), 1],
            top_snp_pph4 = max(my.res$results$SNP.PP.H4),
            status = "ok"
        )

    }, error = function(e) {
        msg <- conditionMessage(e)
        cat(glue("Error processing {gene_name}: {msg}"), "\n")
        return(df_failed(gene_name, msg))
    })
}

# map_df returns a tibble; coerce so the data.table summary below works
df_coloc_res <- as.data.table(map_df(egene_list, run_coloc))

# write results
dir_create(dirname(output_file))
fwrite(df_coloc_res, output_file, row.names = FALSE, sep = "\t")

cat(glue("Saved {nrow(df_coloc_res)} coloc results to {output_file}"), "\n")

# Summarise outcomes so a chromosome producing no usable results is visible in
# the log rather than only as an absent chromosome downstream.
n_ok <- sum(df_coloc_res$status == "ok")
n_fail <- nrow(df_coloc_res) - n_ok
cat(glue("Genes attempted: {length(egene_list)} | succeeded: {n_ok} | failed: {n_fail}"), "\n")
if (n_fail > 0) {
    cat("Failure reasons:\n")
    print(df_coloc_res[status != "ok", .N, by = status][order(-N)])
}
cat("Finished!\n")

close_log()
