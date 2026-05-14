library(fs)
library(tidyverse)
library(qvalue)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

calc_q <- function(p, ...) {
    possibly(qvalue, qvalue_truncp(p))(p, ...) %>% .$qvalues
}

df_cell_type_stats <- map_df(
    map_chr(INPUT, 1) %>% 
        set_names(., str_remove_all(basename(.), "\\..*")),
    ~read_tsv(.x) %>% 
        rename(cell_type = group) %>% 
        mutate(qval_assoc_mcp = calc_q(assoc_mcp),
               bh_assoc_mcp = p.adjust(assoc_mcp, method = "BH")),
    .id = "phenotype"
)

write_tsv(df_cell_type_stats, OUTPUT[[1]])
