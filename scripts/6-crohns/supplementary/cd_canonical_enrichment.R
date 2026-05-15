# enrichment
source("scripts/preprocess.R")

library(data.table)
library(tidyverse)

df_canon <- readRDS("resources/misc/canonical_genes_crohns.rds")

df_crohns <- df_msmr[phenotype == "crohns"]

# 2 x 2 enrichment test
df_test <- data.table(gene = unique(df_msmr_tenk10k$Gene)) %>% 
  mutate(canon = gene %in% df_canon$Gene,
         mr = gene %in% df_crohns$Gene)

table(mr = df_test$mr, canon = df_test$canon) %>% 
  fisher.test(alternative = "greater") %>% 
  broom::tidy()