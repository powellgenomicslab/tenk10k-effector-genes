
# make trait list for TenK10K phase 1

library(readxl)
library(tidyverse)

df_trait_meta <- read_excel("resources/metadata/trait_metadata_curated.xlsx")

df_trait_meta |> 
    filter(include) |> 
    pull(trait_id) |> 
    write_lines(snakemake@output[[1]])