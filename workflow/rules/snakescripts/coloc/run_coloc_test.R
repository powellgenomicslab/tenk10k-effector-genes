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

# Parameters
window_bp <- snakemake@params[["window_bp"]]

# Extract metadata from wildcards
biosample <- snakemake@wildcards[["biosample"]]
pheno <- snakemake@wildcards[["pheno"]]
study <- snakemake@wildcards[["study"]]
chr <- snakemake@wildcards[["chr"]]

# phenotype metadata
df_pheno_meta <- fread(pheno_metadata_file) %>%
    filter(trait_id == pheno)
pheno_type <- ifelse(df_pheno_meta$supercategory == "biological", "quant", "cc")

# assume that quantitative traits have sdY = 1
# todo: add sdY in the metadata file if needed
pheno_sd <- ifelse(pheno_type == "quant", 1, NA)

## 2. Load data -------------------------------------------------------------
df_eqtl <- read_fst(eqtl_file, as.data.table = TRUE)
setkey(df_eqtl, snp)
egene_list <- unique(df_eqtl$gene)

df_gwas <- fread(gwas_file)
setkey(df_gwas, SNP)

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

df_blank <- function(){
    data.table(
        chr = chr,
        biosample = character(0),
        pheno = character(0),
        gene = character(0),
        nsnps_coloc = numeric(0),
        PP.H0.abf = numeric(0),
        PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0),
        PP.H3.abf = numeric(0),
        PP.H4.abf = numeric(0),
        top_snp = character(0),
        top_snp_pph4 = numeric(0)
    )
}

if (nrow(df_gene_loc) == 0) {
    cat(glue("No eGenes found for chromosome {chr} - creating empty output"), "\n")
    # Create empty output with correct structure
    dir_create(dirname(output_file))
    fwrite(df_blank(), output_file, row.names = FALSE, sep = "\t")
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
            top_snp_pph4 = max(my.res$results$SNP.PP.H4)
        )
        
    }, error = function(e) {
        cat(glue("Error processing {gene_name}: {e$message}"), "\n")
        return(df_blank())
    })
}

df_coloc_res <- map_df(egene_list, run_coloc)

# write results
dir_create(dirname(output_file))
fwrite(df_coloc_res, output_file, row.names = FALSE, sep = "\t")

cat(glue("Saved {nrow(df_coloc_res)} coloc results to {output_file}"), "\n")
cat("Finished!\n")
