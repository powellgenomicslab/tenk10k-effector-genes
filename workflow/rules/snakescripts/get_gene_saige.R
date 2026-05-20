library(data.table)
library(tidyverse)
library(fs)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

INPUT <- dir_ls("resources/saige_eqtl/tenk10k_phase1", recurse = TRUE,
                glob = "*common_raw.tsv") %>%
        set_names(., basename(dirname(.))) %>%
        .[names(.) != "CD4_TCM_sample_perm0"] %>% 
        as.list()

get_unique_gene <- function(f, gene_col_n = 18) {
    cmd <- paste0("awk 'NR > 1 {print $", gene_col_n, "}' ", f, " | sort -u")
    fread(cmd = cmd, header = FALSE)
}

df_unique_genes <- map(INPUT, get_unique_gene)

df_wide <- bind_rows(df_unique_genes, .id = "cell_type") %>% 
    mutate(value = TRUE) %>% 
    pivot_wider(names_from = cell_type, values_from = value, values_fill = FALSE) %>%
    mutate(gene = V1)

fwrite(df_wide, "resources/misc/tenk10k_phase1_expression_celltype.tsv", sep = "\t",
       quote = FALSE, row.names = FALSE)

bind_rows(df_unique_genes, .id = "cell_type")  %>% 
    group_by(cell_type) %>%
    tally() %>%
    arrange(desc(n)) %>%
    print.data.frame()