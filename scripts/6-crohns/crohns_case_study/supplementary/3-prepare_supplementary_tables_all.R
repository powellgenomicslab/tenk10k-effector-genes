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

# # data from Bradley's analysis on effector gene enrichment (an intermediate file saved from the TenK10K_causal_inference_vs_IBDverse repo)
# load("resources/crohns_case_study/revision/IBDverse_pi1_tables.RData")

pi1_combined <- read_tsv("results/pi1_mr_novelty_results/pi1_combined.tsv")
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

OR <- read_tsv("results/pi1_mr_novelty_results/tenk_not_novel_res_relative.tsv") 
OR <- OR %>% 
  select(category, nde, n_notnovelde, OR, CI_low, CI_high, pval, prop_blood) %>% 
  rename(`Cell Lineage` = category, 
            `Number of IBDverse Effector Genes` = nde, 
            `Number of Intersecting IBDverse Effector & Crohn's MR Genes` = n_notnovelde, 
            `Odds ratio` = OR,
            `95% lower CI` = CI_low,
            `95% upper CI` = CI_high,
            "Comparison" = pval, 
            "Proportion from blood" = prop_blood)

# Get literature annotation df and change col names for Sup Table
source("scripts/crohns_case_study/supplementary/2-prepare_annot.R")
###########################################################################################

df_list_names <- c("crohns_summary", "crohns_shared_immune", "crohns_mr_deg_comparison", "ibdverse_eqtl_replication", "ibdverse_effectorgene_OR", "crohns_literature_annotations")
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
