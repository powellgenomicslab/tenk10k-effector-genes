library(data.table)
library(arrow)
library(qvalue)
library(tidyr)
library(dplyr)

# Read in trait-level results
input_dir <- "/workspace/ivw-ld"
input_filelist <- Sys.glob(file.path(input_dir, "results/*_IVW-LD_results.tsv"))
input_df_list <- lapply(input_filelist, fread)
output_df <- rbindlist(input_df_list)

# Write to a parquet file for use in the sensitivity pipeline
output_dir <- "/workspace/ivw-ld/results"
output_filepath <- file.path(output_dir, "tenk10k_phase1_sensitivity.ivw-ld.parquet.gz")
write_parquet(output_df, output_filepath, compression = "gzip")

