# extract association results
# Using arrow / dplyr
library(tidyverse)
library(arrow)
library(readxl)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params

# read trait metadata
df_trait_meta <- read_excel(INPUT$trait_metadata)

# helper function to query a dataset

filter_unique <- function(dataset, col, x) {
    x <- discard(unique(x), ~ is.na(.x))
    dataset %>%
        filter({{ col }} %in% x)
}

# Extract target-disease association for available efo_ids
fs::dir_create(PARAMS$output_dir)

assocs <- c("assoc_overall", "assoc_datasource", "assoc_datatype")

for (x in assocs) {
    df_assoc <- list(
        direct = INPUT[[paste0(x, "_direct")]],
        indirect = INPUT[[paste0(x, "_indirect")]]
    ) %>%
        map_df(~open_dataset(.x) %>% 
                filter_unique(diseaseId, df_trait_meta$query_id) %>%
                collect(),
            .id = "association_type"
        )
    write_parquet(df_assoc, OUTPUT[[x]], compression = "gzip")
}