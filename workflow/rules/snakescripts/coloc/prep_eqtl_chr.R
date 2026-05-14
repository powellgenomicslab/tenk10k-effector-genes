# Prepare eQTL data for coloc analysis (per chromosome)
# Converts SAIGE eQTL output to coloc-compatible format

library(tidyverse)
library(data.table)
library(fst)

## Get inputs from snakemake
saige_file <- snakemake@input[["saige"]]
egene_file <- snakemake@input[["egene"]]
output_file <- snakemake@output[[1]]

biosample <- snakemake@wildcards[["biosample"]]
chr <- as.integer(snakemake@wildcards[["chr"]])
study <- snakemake@wildcards[["study"]]

cat(sprintf("Processing eQTL data: %s / %s / chr%d\n", study, biosample, chr))

## Load eGene list
cat("Loading eGene list...\n")
egene_list <- fread(egene_file, header = FALSE)$V1
cat(sprintf("Number of eGenes: %d\n", length(egene_list)))

## Load SAIGE eQTL results
cat("Loading SAIGE eQTL results...\n")
# Assuming columns: variant_id, gene_id, chr, pos, ref, alt, AF, beta, se, pval, etc.
eqtl_df <- fread(saige_file) %>%
    # Filter to current chromosome and eGenes
    filter(chr == !!chr, gene_id %in% egene_list)

if (nrow(eqtl_df) == 0) {
    cat(sprintf("No eQTL data for chromosome %d - creating empty output\n", chr))
    # Create empty data frame with correct structure
    empty_df <- data.table(
        snp = character(0),
        beta = numeric(0),
        varbeta = numeric(0),
        N = numeric(0),
        MAF = numeric(0),
        position = numeric(0),
        gene = character(0)
    )
    write_fst(empty_df, output_file)
    quit(save = "no")
}

cat(sprintf("Found %d eQTL associations\n", nrow(eqtl_df)))

## Convert to coloc format
# Required columns for coloc: snp, beta, varbeta, N, MAF, position
# Optional: gene (for subsetting)
coloc_df <- eqtl_df %>%
    transmute(
        snp = variant_id,  # or paste(chr, pos, ref, alt, sep = ":")
        beta = beta,
        varbeta = se^2,  # variance = standard error squared
        N = n_samples,  # sample size (should be constant per biosample)
        MAF = pmin(AF, 1 - AF),  # minor allele frequency
        position = pos,
        gene = gene_id
    ) %>%
    # Remove any missing values
    filter(!is.na(beta), !is.na(varbeta), !is.na(MAF)) %>%
    as.data.table()

## Save as FST for efficient storage
cat(sprintf("Saving %d associations to %s\n", nrow(coloc_df), output_file))
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_fst(coloc_df, output_file, compress = 50)

cat("Done!\n")
