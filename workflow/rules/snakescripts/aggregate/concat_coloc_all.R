# Aggregate all coloc results into a single parquet file

library(arrow)
library(fs)

## Get inputs from snakemake
input_files <- unlist(snakemake@input)
output_file <- snakemake@output[[1]]
study <- snakemake@wildcards[["study"]]

# Shared schema normalisation - see the header of coloc_schema.R for why the
# per-file column layout varies.
source(snakemake@params[["schema_helper"]])

combined_df <- read_coloc_all_chr(input_files, expected_n = length(input_files))

cat(sprintf("Aggregated %d rows from %d files\n", nrow(combined_df), length(input_files)))

## Save as compressed parquet file
dir_create(dirname(output_file))
write_parquet(combined_df, output_file, compression = "gzip")
