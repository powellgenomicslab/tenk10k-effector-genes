library(tidyverse)

mr_deg <- readRDS("resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS")
colnames(mr_deg)
crohns_annotated_results_deg <- mr_deg %>% 
  #readRDS("resources/crohns_case_study/supplementary/Crohns_Annotated_Results_DEG.RDS") %>% as.data.frame()
  mutate(across(where(is.logical), function(x) replace_na(x, FALSE))) %>% 
  mutate(Category = case_when(
    mr_deg_sign_concord == TRUE & sign(b_SMR) == 1 ~ "MR + DEG +", 
    mr_deg_sign_concord == TRUE & sign(b_SMR) == -1 ~ "MR - DEG -",
    mr_deg_sign_concord == FALSE & sign(b_SMR) == 1 ~ "MR + DEG -",
    mr_deg_sign_concord == FALSE & sign(b_SMR) == -1 ~ "MR - DEG +")) %>% 
  mutate(Concordance = ifelse(mr_deg_sign_concord, "Concordant", "Discordant")) %>% 
  select(Gene, probeID, cell_type, gut_cell_type, major_cell_type, Contrast, Location, Category, Concordance, mr_deg_sign_ct_het, p_transform_mr, p_transform_deg) %>% 
  rename("EnsemblID" = probeID, 
         "Cell Type" = cell_type, 
         "Major Cell Type" = major_cell_type,
         `Mixed Concordance Between Cell Types` = mr_deg_sign_ct_het,
         `Original cell type annotation in Kong, L. et al 2023` = gut_cell_type,
         "log10 P MR x Direction of Effect" = p_transform_mr, 
         "log10 P DEG x Direction of Effect" = p_transform_deg)

#saveRDS(crohns_annotated_results_deg, "resources/crohns_case_study/supplementary/crohns_annotated_results_deg.RDS")