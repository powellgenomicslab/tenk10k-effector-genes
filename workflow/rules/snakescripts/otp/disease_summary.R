# quantify open targets target-disease evidence

# Using arrow / dplyr
library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params

# read trait metadata
# INPUT <- list(
#     otp_disease_dir = "resources/nci/otp_output/25.03/disease"
# )

ds_disease <- open_dataset(INPUT$otp_disease_dir)

df_disease <- ds_disease %>%
    collect()

# get summaries of target-diseae pairs
df_areas <- df_disease %>% 
    filter(ontology$isTherapeuticArea == TRUE) %>% 
    select(therapeuticAreas = id,  therapeuticAreaName = name)

df_disease_summary <- df_disease %>%
    filter(ontology$isTherapeuticArea == FALSE) %>%
    select(id, name, therapeuticAreas) %>% 
    unnest(c(therapeuticAreas)) %>% 
    left_join(df_areas)

fs::dir_create(PARAMS$output_dir)
write_tsv(df_disease_summary, OUTPUT[[1]])
