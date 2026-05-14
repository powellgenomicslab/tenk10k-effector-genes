# Coloc analysis script for Snakemake pipeline
# Modified from original run_coloc.R to work with snakemake

# Load libraries
library(coloc)
library(tidyverse)
library(glue)
library(data.table)
library(fst)

## 1. Get Snakemake parameters -------------------------------------------------------------
# Input files
eqtl_file <- snakemake@input[["eqtl"]]
gwas_file <- snakemake@input[["gwas"]]
egene_file <- snakemake@input[["egene"]]
gene_loc_file <- snakemake@input[["gene_loc"]]

# Output file
output_file <- snakemake@output[["coloc"]]

# Parameters
window_bp <- snakemake@params[["window_bp"]]

# Extract metadata from wildcards
biosample <- snakemake@wildcards[["biosample"]]
pheno <- snakemake@wildcards[["pheno"]]
chr <- snakemake@wildcards[["chr"]]
study <- snakemake@wildcards[["study"]]

cat(glue("Processing: {study} / {biosample} / {pheno} / chromosome {chr}"), "\n")
cat(glue("Window size: {format(window_bp, scientific = FALSE)} bp"), "\n")

## 2. Load data -------------------------------------------------------------
cat("Loading eGene list ...\n")
egene_list <- fread(egene_file, header = FALSE)$V1
cat(glue("Number of eGenes: {length(egene_list)}"), "\n")

cat("Loading gene locations ...\n")
egene_loc_df <- fread(gene_loc_file) %>%
    filter(gene %in% egene_list, chr == paste0("chr", chr)) %>%
    mutate(
        cis_start = start - window_bp,
        cis_end = end + window_bp
    )

if (nrow(egene_loc_df) == 0) {
    cat(glue("No eGenes found for chromosome {chr} - creating empty output"), "\n")
    # Create empty output with correct structure
    empty_result <- data.frame(
        gene = character(0),
        nsnps_coloc_tested = numeric(0),
        PP.H0.abf = numeric(0),
        PP.H1.abf = numeric(0),
        PP.H2.abf = numeric(0),
        PP.H3.abf = numeric(0),
        PP.H4.abf = numeric(0),
        biosample = character(0),
        pheno = character(0),
        chr = character(0),
        top_snp = character(0),
        stringsAsFactors = FALSE
    )
    fwrite(empty_result, output_file, row.names = FALSE)
    quit(save = "no")
}

cat(glue("Number of eGenes for chromosome {chr}: {nrow(egene_loc_df)}"), "\n")

cat("Loading eQTL data ...\n")
eqtl_df <- read_fst(eqtl_file, as.data.table = TRUE)
setkey(eqtl_df, gene)

cat("Loading GWAS data ...\n")
gwas_df <- read_fst(gwas_file, as.data.table = TRUE)
setkey(gwas_df, position)

## 3. Colocalisation analysis -------------------------------------------------------------
result_list <- vector("list", nrow(egene_loc_df))

for (i in seq_len(nrow(egene_loc_df))) {
    gene_name <- egene_loc_df$gene[i]
    
    # Extract cis window coordinates
    cis_start <- egene_loc_df$cis_start[i]
    cis_end <- egene_loc_df$cis_end[i]
    
    # Filter GWAS data based on cis window
    gwas_df_subset <- gwas_df[position %between% c(cis_start, cis_end)]
    
    if (nrow(gwas_df_subset) == 0) {
        cat(glue("No GWAS data for {gene_name} in cis-window: skipping ..."), "\n")
        next
    }
    
    # Process GWAS data for coloc
    # Expected columns: snp, beta, varbeta, N, MAF, position, type, s (for case-control)
    gwas_data <- as.list(gwas_df_subset)
    if (!"type" %in% names(gwas_data)) {
        gwas_data$type <- "cc"  # or "quant" depending on trait
    }
    
    # Filter eQTL data for current gene
    eqtl_data_subset <- eqtl_df[J(gene_name)]
    
    if (nrow(eqtl_data_subset) == 0) {
        cat(glue("No eQTL data for {gene_name}: skipping ..."), "\n")
        next
    }
    
    # Process eQTL data for coloc
    # Expected columns: snp, beta, varbeta, N, MAF, position
    eqtl_data <- as.list(eqtl_data_subset[, !"gene"])
    eqtl_data$type <- "quant"
    
    # Find number of shared SNPs
    shared_snps <- length(intersect(gwas_data$snp, eqtl_data$snp))
    if (shared_snps == 0) {
        cat(glue("No common SNPs between GWAS and eQTL for {gene_name}: skipping ..."), "\n")
        next
    }
    
    cat(glue("Gene {gene_name}: {shared_snps} shared SNPs"), "\n")
    
    # Perform colocalisation analysis
    tryCatch({
        my.res <- coloc.abf(
            dataset1 = eqtl_data,
            dataset2 = gwas_data
        )
        
        # Extract results
        p_df <- data.frame(
            gene = gene_name,
            nsnps_coloc_tested = my.res$summary[1],
            PP.H0.abf = my.res$summary[2],
            PP.H1.abf = my.res$summary[3],
            PP.H2.abf = my.res$summary[4],
            PP.H3.abf = my.res$summary[5],
            PP.H4.abf = my.res$summary[6],
            biosample = biosample,
            pheno = pheno,
            chr = paste0("chr", chr),
            top_snp = arrange(my.res$results, desc(SNP.PP.H4))[1, 1],
            stringsAsFactors = FALSE
        )
        
        result_list[[i]] <- p_df
    }, error = function(e) {
        cat(glue("Error processing {gene_name}: {e$message}"), "\n")
    })
}

## 4. Save results -------------------------------------------------------------
result_df <- rbindlist(result_list, use.names = TRUE, fill = TRUE)

# Create output directory if needed
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Write results
fwrite(result_df, output_file, row.names = FALSE)

cat(glue("Saved {nrow(result_df)} coloc results to {output_file}"), "\n")
cat("Finished!\n")
