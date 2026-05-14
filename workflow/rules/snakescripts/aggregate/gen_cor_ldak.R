library(tidyverse)
library(data.table)

INPUT <- snakemake@input

read_input <- function(input_file) {
  x <- str_split(basename(input_file), "\\.")[[1]]
  fread(input_file) %>%
    mutate(trait1 = x[1], trait2 = x[2], .before = 1)
}

df <- map_df(INPUT, read_input)

fwrite(df, snakemake@output[[1]], sep = "\t", quote = FALSE)