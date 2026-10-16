
# make trait list for TenK10K phase 1

library(readxl)
library(tidyverse)

df_trait_meta <- fread("metadata/trait.tsv")

df_trait_meta |> 
    filter(include) |> 
    pull(trait_id) |> 
    write_lines(snakemake@output[[1]])