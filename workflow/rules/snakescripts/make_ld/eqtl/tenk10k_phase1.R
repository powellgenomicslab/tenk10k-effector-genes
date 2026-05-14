
# Load libraries
library(tidyverse)
library(glue)
library(data.table)
library(fs)
library(fst)
library(sys)

# Interactive
study <- "tenk10k_phase1"
# chr <- 8 
# INPUT <- list(
#     bfile = glue("resources/genotypes/{study}/chr{chr}.{ext}", ext = c("bed", "bim", "fam")),
#     eqtl_dir = glue("resources/coloc/{study}"),
#     gene_loc = glue("resources/misc/gencode.v44.gene_type.tsv")
# )

# OUTPUT <- list(
#         ld = glue("resources/ld/eqtl/{study}/chr{chr}")
# )

# THREADS <- 8
# window_bp <- 1e5
## 1. Get Snakemake parameters -------------------------------------------------------------
# Input files

study <- snakemake@wildcards[["study"]]
chr <- snakemake@wildcards[["chr"]]
INPUT <- snakemake@input
OUTPUT <- snakemake@output
THREADS <- snakemake@threads
window_bp <- snakemake@params[["window_bp"]]

setDTthreads(THREADS)
threads_fst(THREADS)

## 2. Preliminary step --------------------------------------
CHR <- chr
df_gene_loc <- fread(INPUT$gene_loc)  |> 
    filter(chr == CHR) |> 
    mutate(
        cis_start = pmax(1, start - window_bp),
        cis_end = end + window_bp
    )

df_eqtl <- dir_ls(INPUT$eqtl_dir, recurse = TRUE, glob = glue("*chr{chr}.fst")) |> 
    map_dfr(~read_fst(.x, as.data.table = TRUE)) |> 
    distinct(gene, snp) |> 
    mutate(A1 = str_split(snp, ":") |> map_chr(3))

n <- df_gene_loc[,.N]
mk_ld_gene <- function(gene_name){
    bfile_prefix <- glue("resources/genotypes/{study}/chr{chr}")
    cis_start <- df_gene_loc[ensembl_gene_id == gene_name, cis_start]
    cis_end <- df_gene_loc[ensembl_gene_id == gene_name, cis_end]
    df_eqtl_gene <- df_eqtl[gene == gene_name]

    out_prefix <- glue("{OUTPUT$ld}/{gene_name}")

    # write temporary file (for allele matching when creating LD matrix)
    tmp_eqtl <- file_temp()
    fwrite(df_eqtl_gene[,.(A1, snp)], tmp_eqtl, sep = "\t")
     args_ld <- c(
        "--bfile", bfile_prefix,
        "--r", "square",
        "--make-just-bim",
        "--chr", chr,
        "--from-bp", cis_start,
        "--to-bp", cis_end,
        "--a1-allele", tmp_eqtl, "1", "2",
        "--out", out_prefix,
        "--silent"
    )

    i <- df_gene_loc[ensembl_gene_id == gene_name, which = TRUE]

    cat(glue("\nRunning LD matrix creation for {gene_name} ({i}/{n})\n"))

    exec_wait("plink", args_ld)

}

dir_create(OUTPUT$ld)
walk(df_gene_loc$ensembl_gene_id, mk_ld_gene)

cat(glue("Saved {n} LD matrix to {OUTPUT$ld}"), "\n")
cat("Finished!\n")
