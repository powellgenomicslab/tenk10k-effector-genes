###
### This script reads the summary statistics from the TenK10K SMR-multi analysis
### and creates a file of cell types for each trait that has significant results
### (lfdr < 0.05). These files will be used to run IVW-LD and MR-link-2 in the
### sensitivity pipeline
###

args <- commandArgs(trailingOnly = TRUE)
msmr_filepath <- args[1]
output_path <- args[2]
sentinel_file <- args[3]

library(data.table)
library(arrow)
library(qvalue)
library(tidyverse)

# Read in the data and filter out non-significant associations
tenk10k_phase1 <- read_parquet(msmr_filepath, 
    col_select = c("phenotype", "biosample", "probeID", "p_SMR_multi")) |> as.data.table()
tenk10k_phase1[, lfdr_msmr_pheno := qvalue(p_SMR_multi)$lfdr, by = "phenotype"]
tenk10k_phase1 <- tenk10k_phase1[lfdr_msmr_pheno < 0.05]

# Generate list of phenotypes and cell types
tenk10k_phase1 <- tenk10k_phase1[, .(phenotype, biosample)]
tenk10k_phase1 <- unique(tenk10k_phase1)
tenk10k_phase1 <- as_tibble(tenk10k_phase1)

# For each trait, output cell types they are found in
tenk10k_phase1 %>%
  group_by(phenotype) %>%
  group_map(~ write_csv(.x %>% select(biosample), 
    paste0(output_path, .y$phenotype, ".txt"), col_names = FALSE), .keep = FALSE)

# At the end, output a file of phenotypes that have significant results
# Create a file signifying that the script has finished
tenk10k_phenotype_list <- unique(tenk10k_phase1$phenotype)
write_tsv(tibble(phenotype = tenk10k_phenotype_list), sentinel_file, col_names = FALSE)
