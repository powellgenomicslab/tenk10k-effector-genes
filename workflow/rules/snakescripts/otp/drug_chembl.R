# quantify open targets target-disease evidence

# Using arrow / dplyr
library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params

# read trait metadata
# INPUT <- list(
#     otp_evidence_dir = "resources/nci/otp_output/25.06/evidence"
# )

ds <- open_dataset(INPUT$otp_evidence_dir)

df <- ds %>%
    filter(datasourceId == "chembl") %>% 
    select(targetId, diseaseId, score, variantEffect, directionOnTrait) %>%
    collect() %>% 
    distinct()

fs::dir_create(PARAMS$output_dir)
write_tsv(df, OUTPUT[[1]])
