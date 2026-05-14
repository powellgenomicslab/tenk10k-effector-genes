library(arrow)
library(fs)
library(tidyverse)

files <- snakemake@input

read_data <- function(x) {
    pheno <- basename(x) %>% str_remove_all("\\.magma\\.tsv")
    read_tsv_arrow(x) %>% 
        mutate(phenotype = pheno,
               .before = GENE)
}

map_df(files, read_data) %>% 
    write_parquet(
        snakemake@output[[1]],
        compression = "gzip"
    )
