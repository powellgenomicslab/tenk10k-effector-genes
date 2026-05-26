args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_filepath <- args[2]

library(data.table)
library(dplyr)
library(tidyr)
library(arrow)
library(fs)

setDTthreads(4)

# Grab list of aggregated trait output files
input_filelist <- Sys.glob(file.path(input_dir, "*/*_mrlink2_aggregated.tsv"))

# Prepare list of files as a data.table, to split the strings
metadata_dt <- data.table(filepath = input_filelist)
metadata_dt[, cell_type := basename(dirname(filepath))]
metadata_dt[, trait := sub("_mrlink2_aggregated\\.tsv$", "", basename(filepath))]

input_dt_list <- lapply(metadata_dt[["filepath"]], fread)
names(input_dt_list) <- metadata_dt[["filepath"]]
input_dt <- rbindlist(input_dt_list, idcol = "filepath")
output_dt <- merge(input_dt, metadata_dt, by = "filepath", all.x = TRUE)
output_dt[, filepath := NULL]

write_parquet(output_dt, output_filepath, compression = "gzip")