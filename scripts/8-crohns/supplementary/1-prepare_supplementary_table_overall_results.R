# Supplementary Table
library(tidyverse)
library(writexl)
library(formattable)

mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS")

all_res <- mr %>% group_by(cell_type) %>% tally() %>% ungroup() %>% rename(All = n)

condition <- mr %>% group_by(eqtlgen_mr, cell_type) %>% tally() %>% filter(eqtlgen_mr == F) %>% ungroup() %>% select(-eqtlgen_mr) %>% rename(`Not in Bulk MR` = n)
both <- left_join(all_res, condition, by = "cell_type") %>% 
  mutate(`Prop Not in BulkMR` = `Not in Bulk MR`/All)

bulk_eqtl_comparison <- both


condition <- mr %>% group_by(magma_gene, cell_type) %>% tally() %>% filter(magma_gene == F) %>% ungroup() %>% select(-magma_gene) %>% rename(`Not in MAGMA` = n)
both <- left_join(all_res, condition, by = "cell_type") %>% 
  mutate(`Prop Not in MAGMA` = `Not in MAGMA`/All)

magma_comparison <- both

#########################################################################################################

condition <- mr %>% group_by(ct_spec, cell_type) %>% tally() %>% filter(ct_spec == T) %>% ungroup() %>% select(-ct_spec) %>% rename(`Cell Type Specific` = n)
both <- left_join(all_res, condition, by = "cell_type") %>% 
  mutate(`Prop Cell Type Specific` = `Cell Type Specific`/All)

ct_spec_comparison <- both

# Combine all 3 dfs 
df_list <- list(bulk_eqtl_comparison, magma_comparison, ct_spec_comparison)

combined <- df_list %>%
  reduce(full_join, by = c("cell_type", "All"))


# Prop cell type specific not in eQTL Gene
mr_spec <- mr %>% filter(ct_spec == T) %>% group_by(cell_type) %>% tally() %>% rename(`Cell Type Specific` = n)
condition <- mr %>% filter(ct_spec == T) %>% group_by(eqtlgen_mr, cell_type) %>% tally() %>% filter(eqtlgen_mr == F) %>% ungroup() %>% select(-eqtlgen_mr) %>% rename(`Cell Type Specific Not in Bulk` = n)
both <- left_join(mr_spec, condition, by = "cell_type") %>% 
  mutate(`Prop Not in Bulk (Of Cell Type Specific)` = `Cell Type Specific Not in Bulk`/`Cell Type Specific`)

ct_spec_bulk_comparison <- both

# Combine with other dataframe, note the duplicated Cell Type Specific Column
combined <- full_join(combined, ct_spec_bulk_comparison, by = c("cell_type", "Cell Type Specific")) %>% 
  mutate(`Prop Not in Bulk (Of Non-Cell Type Specific)` = (`Not in Bulk MR`-`Cell Type Specific Not in Bulk`)/(All-`Cell Type Specific`))

# get the mean of every column exept the first 
combined_summary <- combined %>%
  ungroup() %>% 
  summarise(across(-1, mean, na.rm = TRUE))

index <- grep("Prop", colnames(combined))
index2 <- grep("Prop", colnames(combined_summary))

# replace NAs with 0 and make the proportions percentages 
combined <- combined %>%
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  mutate(across(index, ~ percent(., digits = 1)))

combined_summary <- combined_summary %>%
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  mutate(across(index2, ~percent(., digits = 1)))

# name the tables as sheet names 
tables <- list(Crohns_Summary = combined, Crohns_Summary_Mean = combined_summary)

# note that the percentage format is not preserved in excel which is annoying. I tried and tested this. The file was changed after export.  
# writexl::write_xlsx(tables, "resources/crohns_case_study/supplementary/S11_Crohns_Summary.xlsx")

save(combined, combined_summary, file = "resources/crohns_case_study/supplementary/crohns_summary_gene_count_df.RData")








