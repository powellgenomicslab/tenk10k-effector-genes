## Purpose: Summarize Open Targets drug evidence from ChEMBL
## Input: Open Targets evidence dataset directory
## Output: Drug evidence summary TSV file

library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params


ds <- open_dataset(INPUT$otp_evidence_dir)

df <- ds %>%
    filter(datasourceId == "chembl") %>% 
    select(targetId, diseaseId, score, variantEffect, directionOnTrait) %>%
    collect() %>% 
    distinct()

fs::dir_create(PARAMS$output_dir)
write_tsv(df, OUTPUT[[1]])
