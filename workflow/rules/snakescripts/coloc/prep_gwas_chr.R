# Prepare GWAS data for coloc analysis (per chromosome)
# Converts .ma format GWAS summary stats to coloc-compatible format

library(tidyverse)
library(data.table)
library(fst)

## Get inputs from snakemake
ma_file <- snakemake@input[["ma"]]
output_file <- snakemake@output[[1]]

pheno <- snakemake@wildcards[["pheno"]]
chr <- as.integer(snakemake@wildcards[["chr"]])

cat(sprintf("Processing GWAS data: %s / chr%d\n", pheno, chr))

## Load GWAS summary statistics
cat("Loading GWAS summary stats...\n")
# .ma format columns: SNP, A1, A2, freq, b, se, p, n
gwas_df <- fread(ma_file)

## Extract chromosome and position from SNP column
# Assuming SNP format is chr:pos:ref:alt or similar
cat("Extracting chromosome and position...\n")
gwas_df <- gwas_df %>%
    separate(SNP, into = c("snp_chr", "snp_pos", "snp_ref", "snp_alt"), 
             sep = ":", remove = FALSE, convert = TRUE) %>%
    filter(snp_chr == chr)

if (nrow(gwas_df) == 0) {
    cat(sprintf("No GWAS data for chromosome %d - creating empty output\n", chr))
    # Create empty data frame with correct structure
    empty_df <- data.table(
        snp = character(0),
        beta = numeric(0),
        varbeta = numeric(0),
        N = numeric(0),
        MAF = numeric(0),
        position = numeric(0),
        type = character(0),
        s = numeric(0)
    )
    write_fst(empty_df, output_file)
    quit(save = "no")
}

cat(sprintf("Found %d SNPs\n", nrow(gwas_df)))

## Convert to coloc format
# Required columns for coloc: snp, beta, varbeta, N, MAF, position, type
# For case-control: also need 's' (proportion of cases)
coloc_df <- gwas_df %>%
    transmute(
        snp = SNP,
        beta = b,
        varbeta = se^2,  # variance = standard error squared
        N = n,  # sample size
        MAF = pmin(freq, 1 - freq),  # minor allele frequency
        position = snp_pos,
        type = "cc"  # "cc" for case-control, "quant" for quantitative
        # For case-control traits, uncomment and set appropriately:
        # s = n_cases / n  # proportion of cases
    ) %>%
    # Remove any missing values
    filter(!is.na(beta), !is.na(varbeta), !is.na(MAF), !is.na(N)) %>%
    as.data.table()

# Note: You may need to add trait-specific metadata to determine:
# 1. Whether trait is case-control or quantitative (type)
# 2. For case-control, the proportion of cases (s)

## Save as FST for efficient storage
cat(sprintf("Saving %d SNPs to %s\n", nrow(coloc_df), output_file))
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_fst(coloc_df, output_file, compress = 50)

cat("Done!\n")
