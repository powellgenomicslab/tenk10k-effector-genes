library(tidyverse)
library(writexl)

load("resources/crohns_case_study/supplementary/crohns_summary_gene_count_df.RData")

crohns_summary <- combined %>% 
  as.data.frame() %>% 
  mutate(across(where(is.numeric), function(x) round(x, digits = 2)))

rm(combined)
rm(combined_summary)

#########################################################################################################################
source("scripts/crohns_case_study/supplementary/2-prepare_crohns_annotated_results.R")
source("scripts/crohns_case_study/supplementary/2-prepare_supplementary_table_cd_genes_in_immune_diseases.R")
source("scripts/crohns_case_study/supplementary/2-prepare_supplementary_table_mr_deg.R")
load("resources/crohns_case_study/revision/IBDverse_pi1_tables.RData")
# pi1_results <- pi1_summary
# data used to generate fig 6e
library(data.table)
pi1_results <- pi1_combined %>% 
  select(discovery, pi1, label_new, major_cell_type, comparison) %>% 
  left_join(pi1_summary %>% select(-max_pi1, -max_dataset_id), by = "discovery") %>% 
  select(discovery, label_new, pi1, min_pi1, median_pi1, comparison, major_cell_type) %>% 
  rename(`TenK10K Cell Type` = discovery, 
         `Max pi1 IBDverse Cell Type` = label_new, 
         `TenK10K Major Cell Type` = major_cell_type, 
         `Max pi1` = pi1,
         `Min pi1` = min_pi1,
         `Median pi1` = median_pi1,
         "Comparison" = comparison)

# data used to generate sup figure for pi1 analysis
# sup_figure_pi1 <- pi1_all[gut == TRUE] %>% 
#   select(discovery, neqtls_intersection, pi1, tissue_label, label_new)

# data for figure 6e
OR <- readRDS("resources/crohns_case_study/revision/OR_IBDverse_coloc_genes_in_CD_MR_genes.RDS")
source("scripts/crohns_case_study/supplementary/2-prepare_annot.R")
###########################################################################################

df_list_names <- c("15-Annotated CD MR Associations", "16-CD Associations in Immune", "17-CD MR-DEG Comparison", "18-IBDverse eQTL Replication", "19-IBDverse-MR Effector Gene OR", "20-All Literature Annotations")
#df_list <- list(crohns_summary, crohns_annotated_results, crohns_immune, crohns_annotated_results_deg, pi1_results, annot)
df_list <- list(crohns_annotated_results, crohns_immune, crohns_annotated_results_deg, pi1_results, OR, annot)

names(df_list) <- df_list_names

names(df_list)

writexl::write_xlsx(df_list, "resources/crohns_case_study/supplementary/Crohns_Summary.xlsx")

df_colnames_list <- map(
  df_list,
  ~tibble(name = names(.x), label = names(.x), description =NA)
)

writexl::write_xlsx(df_colnames_list, "resources/crohns_case_study/supplementary/Crohns_Summary_colnames.xlsx")
