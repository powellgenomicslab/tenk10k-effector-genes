# Extract gene positions from GTF file for coloc analysis
# This creates the gene_pos.csv file required by the coloc pipeline

library(tidyverse)
library(data.table)
library(rtracklayer)

## Get inputs from snakemake
gtf_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
study <- snakemake@wildcards[["study"]]

cat(sprintf("Extracting gene positio11ns for study: %s\n", study))
cat(sprintf("GTF file: %s\n", gtf_file))

## Read GTF file
cat("Reading GTF file...\n")
gtf <- rtracklayer::import(gtf_file)

## Convert to data frame and filter to genes
cat("Extracting gene annotations...\n")
genes <- as.data.frame(gtf) %>%
    filter(type == "gene") %>%
    select(gene_id, seqnames, start, end, gene_name, gene_type) %>%
    rename(
        gene = gene_id,
        chr = seqnames
    ) %>%
    # Ensure chromosome format is consistent (chr1, chr2, etc.)
    mutate(chr = ifelse(str_starts(chr, "chr"), chr, paste0("chr", chr))) %>%
    # Filter to main chromosomes
    filter(chr %in% paste0("chr", c(1:22, "X", "Y"))) %>%
    arrange(chr, start) %>%
    as.data.table()

cat(sprintf("Extracted %d genes\n", nrow(genes)))

## Summary by chromosome
summary_by_chr <- genes %>%
    group_by(chr) %>%
    summarise(n_genes = n()) %>%
    arrange(chr)

cat("\nGenes per chromosome:\n")
print(summary_by_chr)

## Summary by gene type
if ("gene_type" %in% names(genes)) {
    summary_by_type <- genes %>%
        group_by(gene_type) %>%
        summarise(n_genes = n()) %>%
        arrange(desc(n_genes)) %>%
        head(10)

    cat("\nTop 10 gene types:\n")
    print(summary_by_type)
}

## Save gene positions
cat(sprintf("\nSaving gene positions to %s\n", output_file))
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Save with essential columns
gene_pos <- genes %>%
    select(gene, chr, start, end) %>%
    arrange(chr, start)

fwrite(gene_pos, output_file)

cat(sprintf("Saved %d gene positions\n", nrow(gene_pos)))
cat("Done!\n")
