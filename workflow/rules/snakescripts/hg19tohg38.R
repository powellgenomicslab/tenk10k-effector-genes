
library(data.table)
library(tidyverse)
library(rtracklayer)
library(GenomeInfoDb)
library(GenomicRanges)

# Load HRC data
liftover2hg38 <- function(x, chain_file){
    chain <- rtracklayer::import.chain(chain_file)
    
    # Liftover to hg38 position
    # x should be a dataframe formatted in the following format:
    # data.table(CHROM = chromosome, POS = position (hg19), snpid = SNP identifier to link back to the original gwas)
    x <- as.data.table(x)
    x <- x[, end := start]
    x <- x[, strand := "*"]

    # Reorder columns
    x <- x[, .(seqnames, start, end, strand, snpid)]
    coord_df <- as.data.frame(x)
    coord_granges <- GRanges(coord_df)
    genome(coord_granges) <- "hg19"

    # Liftover to hg19
    seqlevelsStyle(coord_granges) = "UCSC"  # necessary
    coord_granges_hg38 = rtracklayer::liftOver(coord_granges, chain)
    coord_granges_hg38 = unlist(coord_granges_hg38)
    genome(coord_granges_hg38) = "hg38"

    # Back to a data table
    hg38_coords <- as.data.table(data.frame(coord_granges_hg38))
    hg38_coords <- hg38_coords[, seqnames := gsub("chr", "", seqnames)]
    hg38_coords <- hg38_coords[, .(seqnames, start, snpid)]
    colnames(hg38_coords) <- c("chr", "pos_hg38", "snpid")
    hg38_coords[, chr := as.numeric(chr)]
    return(hg38_coords)
}