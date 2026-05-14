# Aggregate all coloc results into a single parquet file

library(tidyverse)
library(data.table)
library(arrow)
library(fs)

## Get inputs from snakemake
# This is a list of all coloc result files
input_files <- snakemake@input
output_file <- snakemake@output[[1]]
study <- snakemake@wildcards[["study"]]

# read_coloc_input <- function(file) {
#   df <- fread(file)
#   if ("nspns_coloc_tested" %in% colnames(df)) {
#     df <- df %>% rename(nspns_coloc = nsnps_coloc_tested)
#   }
#   return(df)
# }

# input_files <- dir_ls("results/coloc/tenk10k_phase1", recurse = TRUE, glob = "*all_chr.coloc.tsv")
# crohns_files <- input_files[str_detect(input_files, "crohns")]
# input_files <- dir_ls("results/coloc/tenk10k_phase1", recurse = TRUE, glob = "*all_chr.coloc.tsv")

combined_df <- rbindlist(map(input_files, ~fread(.x, fill = TRUE)), fill = TRUE)
combined_df[is.na(nsnps_coloc_tested), nsnps_coloc_tested := nsnps_coloc]
# combined_df[is.na(nsnps_coloc_tested), nsnps_coloc_tested := V13]

# remove redundant column
combined_df[,nsnps_coloc := NULL]
# combined_df[,V13 := NULL]
## Summary statistics

## Save as compressed parquet file
# output_file <- "results/aggregate/tenk10k_phase1.coloc.parquet.gz"
dir_create(dirname(output_file))
write_parquet(combined_df, output_file, compression = "gzip")
