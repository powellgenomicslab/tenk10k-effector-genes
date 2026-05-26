args <- commandArgs(trailingOnly = TRUE)
#args <- c("t2dm")
trait <- args[1]

library(data.table)
setDTthreads(7)

input_dir <- "/workspace/ivw-ld/output"
input_filelist <- Sys.glob(file.path(input_dir, "*", paste0(trait, "*_IVW-LD_results.tsv")))

input_df_list <- lapply(input_filelist, fread)
output_df <- rbindlist(input_df_list)
output_filepath <- file.path("/workspace/ivw-ld/results", paste0(trait, "_IVW-LD_results.tsv"))
fwrite(output_df, output_filepath, sep = "\t")
