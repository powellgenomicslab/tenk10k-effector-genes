library(tidyverse)
source("scripts/util/write_table.R")

load("resources/crohns_case_study/supplementary/crohns_summary_gene_count_df.RData")

crohns_summary <- combined %>% 
  as.data.frame() %>% 
  mutate(across(where(is.numeric), function(x) round(x, digits = 2)))

rm(combined)
rm(combined_summary)

#########################################################################################################################
source("scripts/8-crohns/supplementary/2-prepare_crohns_annotated_results.R")
source("scripts/8-crohns/supplementary/2-prepare_supplementary_table_cd_genes_in_immune_diseases.R")
source("scripts/8-crohns/supplementary/2-prepare_supplementary_table_mr_deg.R")

# data from Bradley's analysis on effector gene enrichment (an intermediate file saved from the TenK10K_causal_inference_vs_IBDverse repo)
load("resources/crohns_case_study/revision/IBDverse_pi1_tables.RData")
# pi1_results <- pi1_summary
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

# data from Bradley's analysis on effector gene enrichment (an intermediate file saved from the TenK10K_causal_inference_vs_IBDverse repo)
OR <- readRDS("resources/crohns_case_study/revision/OR_IBDverse_coloc_genes_in_CD_MR_genes.RDS") %>% 
  select(category, nde, n_notnovelde, OR, CI_low, CI_high, pval, prop_blood)
# column labels (pval -> "P", etc.) come from the ibdverse_effectorgene_OR sheet of
# metadata/table_column_names.xlsx via write_table()

source("scripts/8-crohns/supplementary/2-prepare_annot.R")
###########################################################################################

df_list_names <- c("crohns_summary", "crohns_shared_immune", "crohns_mr_deg_comparison", "ibdverse_eqtl_replication", "ibdverse_effectorgene_OR", "crohns_literature_annotations")
#df_list <- list(crohns_summary, crohns_annotated_results, crohns_immune, crohns_annotated_results_deg, pi1_results, annot)
df_list <- list(crohns_annotated_results, crohns_immune, crohns_annotated_results_deg, pi1_results, OR, annot)

names(df_list) <- df_list_names

names(df_list)

# Supplementary Tables 16-21, in the order of df_list_names; column labels come from the
# sheet of the same name in metadata/table_column_names.xlsx
for (i in seq_along(df_list)) {
  write_table(df_list[[i]], df_list_names[i], 15 + i)
}
