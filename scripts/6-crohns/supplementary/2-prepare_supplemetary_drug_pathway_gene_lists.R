load("resources/crohns_case_study/crohns_relevant_gene_lists/prepare_drug_pathway_genes.RData")

df1 <- as.data.frame(il12_il23)
df1$pathway_id <- "PID_IL12_2PATHWAY_PID_IL23_PATHWAY"
colnames(df1) <- c("Gene", "Pathway_ID_Combinations")

df2 <- as.data.frame(integrin)
df2$pathway_id <- "PID_INTEGRIN5_PATHWAY_PID_INTEGRIN1_PATHWAY_PID_INTEGRIN_CS_PATHWAY"
colnames(df2) <- c("Gene", "Pathway_ID_Combinations")

df3 <- as.data.frame(tnf)
df3$pathway_id <- "PID_TNF_PATHWAY"
colnames(df3) <- c("Gene", "Pathway_ID_Combinations")

df4 <- as.data.frame(drug_target_genes)
df4$pathway_id <- "NA_Drug_Targets_from_Literature"
colnames(df4) <- c("Gene", "Pathway_ID_Combinations")

all <- rbind(df1, df2, df3, df4)

library(tidyverse)
all_wider <- all %>% 
  mutate(value = TRUE) %>% 
  pivot_wider(id_cols = Gene, names_from = "Pathway_ID_Combinations", values_fill = FALSE) 

write.csv(all_wider, "resources/crohns_case_study/supplementary/drug_pathway_genes_df.csv", row.names = F)
saveRDS(all_wider, "resources/crohns_case_study/supplementary/drug_pathway_genes_df.RDS")