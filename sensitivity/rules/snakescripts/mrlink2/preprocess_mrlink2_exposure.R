###
### Extract significant cell type, chromosome and genes from SMR results regardless of 
### traits from SAIGE-QTL outputs
### 

args <- commandArgs(trailingOnly = TRUE)
msmr_filepath <- args[1]
cell_type <- args[2]
output_dir <- args[3]

library(data.table)
library(arrow)
library(dplyr)
library(fs)
library(qvalue)

setDTthreads(4)

writeExposure <- function(gene, chr_output_dir = NULL, output_dt = NULL){
    gene_dt <- output_dt[probeID == gene]
    gene_dt[, probeID := NULL]
    output_filepath <- file.path(chr_output_dir, paste0(gene, ".txt"))
    fwrite(gene_dt, output_filepath, sep = "\t")
}

processChrom <- function(chrom, output_dt = NULL, output_dir = NULL){
    # Create directory for chromosome if it doesn't exist
    chr_output_dir <- file.path(output_dir, paste0("chr", chrom))
    chr_output_dt <- output_dt[chromosome == chrom]
    dir_create(chr_output_dir, recurse = TRUE)
    
    # For each chromosome, iterate over list of genes and return list of genes and their regions to pass on to Snakemake
    chrom_gene_list <- unique(chr_output_dt$probeID)

    for (gene in chrom_gene_list){
        writeExposure(gene, chr_output_dir = chr_output_dir, output_dt = chr_output_dt)
    }
}


# Load MSMR results and filter for significant associations for that cell type
dt <- read_parquet(msmr_filepath, as_data_frame = FALSE) |> 
    filter(lfdr_msmr_pheno < 0.05, biosample == cell_type) |> as.data.table()

# Get cell types, chromosomes, and gene IDs and make it distinct
dt <- dt[, .(chromosome, probeID)]
dt <- unique(dt)

# For each cell type and chromosome, read in their corresponding eQTL data 
# and prepare exposure files
chrom_list <- unique(dt$chromosome)

# Read in eQTL data
eqtl_filepath <- file.path("/eqtl", cell_type, paste0(cell_type, "_common_all_cis_raw_pvalues.tsv"))
eqtl_dt <- fread(eqtl_filepath)

# Rename eQTL columns to match merge keys
setnames(eqtl_dt, c("MarkerID", "CHR", "gene"), c("variant", "ProbeChr", "probeID"))

# Filter to genes significant in mSMR (right join on chromosome + gene)
probe_dt <- copy(dt)
setnames(probe_dt, "chromosome", "ProbeChr")
eqtl_dt <- merge(eqtl_dt, probe_dt, by = c("ProbeChr", "probeID"), all.y = TRUE)

# Calculate z-score and filter for genome-wide significant eQTL variants
eqtl_dt[, z := BETA / SE]
eqtl_dt <- eqtl_dt[`p.value` < 5e-8]

# Select and rename to GWAS catalogue format
# rsid,chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,z,p_value,n
output_dt <- eqtl_dt[, .(probeID, variant, ProbeChr, POS, Allele2, Allele1, BETA, SE, z, `p.value`, N)]
setnames(output_dt,
    c("variant", "ProbeChr", "POS", "Allele2", "Allele1", "BETA", "SE", "z", "p.value", "N"),
    c("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "z", "p_value", "n"))

for (chrom in chrom_list){
    processChrom(chrom, output_dt = output_dt, output_dir = output_dir)
}
