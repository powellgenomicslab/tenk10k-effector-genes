args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
input_filepath <- args[2]
output_filepath <- args[3]

library(data.table)
library(arrow)
library(dplyr)

setDTthreads(7)

# GWAS format: rsid,chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,z,p_value,n
# ma format: SNP,A1,A2,freq,b,se,p,N

# Read in MA format
gwas_df <- fread(input_filepath)

# Retrieve chromosome and BP location from SNP
gwas_df[, c("chromosome", "base_pair_location") := tstrsplit(SNP, ":", keep = c(1, 2))]

# Calculate z-score
gwas_df[, z := b/se]

# Reorder columns to match GWAS format
# Catch all for N
if ("N" %in% colnames(gwas_df)) {
    setnames(gwas_df, "N", "n")
}

gwas_df <- gwas_df[, c("SNP", "chromosome", "base_pair_location", "A1", "A2", "b", "se", "z", "P", "n")]
output_df <- gwas_df[, .(rsid = SNP, chromosome, base_pair_location, effect_allele = A1, other_allele = A2, beta = b, standard_error = se, z, p_value = P, n = n)]
output_df[, c("chromosome", "base_pair_location") := list(as.integer(chromosome), as.integer(base_pair_location))]

# Save as parquet
write_parquet(output_df, output_filepath, compression = "gzip")