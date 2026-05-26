args <- commandArgs(trailingOnly = TRUE)
#args <- c("lymph", "NK")

trait <- args[1]
cell_type <- args[2]
msmr_filepath <- args[3]
output_filepath <- args[4]

library(data.table)
library(arrow)
library(dplyr)
library(fs)
library(qvalue)

setDTthreads(7)

dt <- read_parquet(msmr_filepath, as_data_frame = FALSE) |>
    filter(phenotype == trait, biosample == cell_type, lfdr_msmr_pheno < 0.05) |> as.data.table()

# For each trait and cell type, produce a file with chromosomes and genes to test in MR-link-2
dt <- dt[, .(chromosome, probeID)]
dt <- unique(dt)

# Write a regions file
fwrite(dt, output_filepath, sep = "\t")
