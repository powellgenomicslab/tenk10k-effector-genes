# quantify open targets target-disease evidence

# Using arrow / dplyr
library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params

# read trait metadata
# INPUT <- list(
#     otp_drug_dir = "resources/nci/otp_output/25.03/known_drug"
# )

ds_drug <- open_dataset(INPUT$otp_drug_dir)

df_drug <- ds_drug %>%
    collect()

# get summaries of target-diseae pairs
df_drug_summary <- df_drug %>%
    mutate(n_trial = map_dbl(urls, length)) %>%
    select(where(~ !is.list(.x)))

fs::dir_create(PARAMS$output_dir)
write_tsv(df_drug_summary, OUTPUT$drug_summary)
