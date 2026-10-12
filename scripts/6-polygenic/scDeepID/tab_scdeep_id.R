# add scdeepID results

source("scripts/util/write_table.R")

library(data.table)
library(tidyverse)

# the pathway label was manually created
# df_res <- fread("results/scDeepID/pathway_test.csv") %>% 
#   mutate(pathway_label = str_remove_all(pathway, "^GOBP_") %>% 
#            str_replace_all("_", " ") %>%
#            str_to_sentence())

# df_pathway_label <- fread("resources/misc/gobp_pathway_label.tsv")
df_pathway_label <- fread("resources/enrichment/gmt/GO_bp_2024_terms.tsv") |> 
  rename(pathway = term_id, pathway_label = term_name)

df_res <- fread("results/scDeepID/pathway_test_full.csv") %>%
  left_join(df_pathway_label, by = c("pathway" = "pathway"))

write_table(df_res, "scdrs_pathway", 10)
