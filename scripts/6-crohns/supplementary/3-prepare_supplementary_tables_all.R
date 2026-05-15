library(tidyverse)
library(writexl)

load("resources/crohns_case_study/supplementary/crohns_summary_gene_count_df.RData")

crohns_summary <- combined %>% as.data.frame()
#crohns_summary_mean <- combined_summary %>% as.data.frame()
crohns_annotated_results <- readRDS("resources/crohns_case_study/supplementary/Crohns_Annotated_Results.RDS") %>% as.data.frame()
crohns_annotated_results_deg <- readRDS("resources/crohns_case_study/supplementary/Crohns_Annotated_Results_DEG.RDS") %>% as.data.frame()

canonical_genes <- read.csv("resources/crohns_case_study/supplementary/canonical_gene_info_table_etal_2023.csv") %>% as.data.frame()
drug_target_or_pathway_genes <- readRDS("resources/crohns_case_study/supplementary/drug_pathway_genes_df.RDS") %>% as.data.frame()
#harmonised_annotation_MR_DEG <- read.csv("resources/crohns_case_study/deg/cell_annotation_or_features/deg_celltype_groups.csv") %>% as.data.frame()
harmonised_annotation_MR_DEG_h5ad <- read.csv("resources/crohns_case_study/deg/cell_annotation_or_features/adata_colon_celltype_groups.csv") %>% as.data.frame()

df_list_names <- c("crohns_summary", "crohns_annotated_results","crohns_annotated_results_deg", "canonical_genes", "drug_target_or_pathway_genes", "harmonised_annot_MR_DEG_h5ad")

df_list <- list(crohns_summary, crohns_annotated_results, crohns_annotated_results_deg, canonical_genes, drug_target_or_pathway_genes, harmonised_annotation_MR_DEG_h5ad)

names(df_list) <- df_list_names

names(df_list)

writexl::write_xlsx(df_list, "resources/crohns_case_study/supplementary/Crohns_Summary.xlsx")

df_colnames_list <- map(
  df_list,
  ~tibble(name = names(.x), label = names(.x), description =NA)
)

writexl::write_xlsx(df_colnames_list, "resources/crohns_case_study/supplementary/Crohns_Summary_colnames.xlsx")

