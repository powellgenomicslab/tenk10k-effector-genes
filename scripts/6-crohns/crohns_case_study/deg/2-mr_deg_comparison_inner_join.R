# Purpose: Intersect DEG and MR genes 
# Take the smallest p-value for a given major cell type - cell type pair 
library(here)
library(tidyverse) 

deg <- readRDS("resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS") %>% 
  select(Gene, Contrast, gut_cell_type, everything()) %>% 
  filter(`Discrete FDR` < 0.05) %>% 
  mutate(deg_sig = TRUE) %>% 
  group_by(Gene, major_cell_type) %>% 
  slice_min(`Discrete DE coefficients p value`, n = 1) %>% 
  ungroup() %>% 
  mutate(p_transform = -log10(`Discrete DE coefficients p value`)*sign(`Discrete DE coefficients`)) 

# update MR major cell type category 
mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% 
  #filter(mr_sens_coloc) %>% 
  mutate(major_cell_type = as.character(major_cell_type) %>% dplyr::replace_when(cell_type == "ILC" ~ "NK", cell_type == "Plasmablast" ~ "Plasma B")) %>% 
  group_by(Gene, major_cell_type) %>% 
  slice_min(p_SMR_multi, n = 1) %>% 
  ungroup() %>% 
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR))

print("Number of significant DEGs: ")
length(unique(deg$Gene))

# filter deg genes
# deg_filtered <- deg %>% filter(Gene %in% mr$Gene) 
# length(unique(deg_filtered$Gene))
# 135 
print("Number of intersecting MR + sig DEGs (probeID and Gene: ")
deg %>% filter(probeID %in% mr$probeID) %>% pull(probeID) %>% unique() %>% length() 
deg %>% filter(Gene %in% mr$Gene) %>% pull(Gene) %>% unique() %>% length() 

# Extra probe IDs, probs non-coding 
############################################################################################################################################################
library(ggrepel) 
print("Inner join MR and DEG results by gene and major cell type")

# Merged the MR and DEG dataframes on gene and major cell type
mr_deg <- inner_join(mr, deg, by = c("Gene", "probeID", "major_cell_type"), suffix = c("_mr", "_deg")) %>% 
  mutate(deg_maj_ct = TRUE) %>% 
  mutate(mr_deg_sign_concord = sign(b_SMR) == sign(`Discrete DE coefficients`)) %>% 
  group_by(Gene) %>% 
  mutate(mr_deg_sign_ct_het = ifelse(n_distinct(mr_deg_sign_concord) == 1, FALSE, TRUE)) %>% 
  ungroup() %>% 
  select(Gene, probeID, major_cell_type, cell_type, gut_cell_type, deg_sig, deg_maj_ct, mr_deg_sign_concord, mr_deg_sign_ct_het, everything()) 

print(mr_deg)

# Print statistics for writing 

print("Number of unique genes MR significant results: ")
length(unique(mr$Gene))

print("Number of unique genes deg filtered to significant results: ")
length(unique(deg$Gene))

print("Number of rows in mr_deg: ")
n_distinct(mr_deg) 

print("Number of unique genes in mr_deg: ")
length(unique(mr_deg$Gene))

print("Number of concordant gene-pbmc cell type-gut cell type combos: ")
mr_deg %>% group_by(mr_deg_sign_concord) %>% tally()

print("Among the concordant mr_deg results at the major cell type level, how many were discordant between cell types?")
mr_deg %>% filter(mr_deg_sign_concord == TRUE) %>% group_by(mr_deg_sign_ct_het) %>% tally()

# how many genes have concordant directions of effect (true or false), how many genes have both concordant and discordant effects depending on cell types?
sum_df <- mr_deg %>% 
  group_by(Gene) %>% 
  summarize(all_true = all(mr_deg_sign_concord), all_false = all(!mr_deg_sign_concord)) %>% 
  group_by(all_true, all_false) %>% 
  tally() 

sum_df

# save to deg
saveRDS(mr_deg, "resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS")

# save text files for pathway enrichment 
mr_deg_con <- mr_deg %>% filter(mr_deg_sign_concord == TRUE) %>% pull(Gene) %>% unique()
mr_deg_discon <- mr_deg %>% filter(mr_deg_sign_concord == FALSE) %>% pull(Gene) %>% unique()
mr_deg_mixed <- mr_deg %>% filter(mr_deg_sign_ct_het == TRUE) %>% pull(Gene) %>% unique()

gene_list <- list(mr_deg_con = mr_deg_con, mr_deg_discon = mr_deg_discon, mr_deg_mixed = mr_deg_mixed) 

print("Print the gene lists for concordant, discordant, mixed: ")
gene_list

for (i in seq_along(gene_list)){
  df <- gene_list[[i]]
  write.table(df, file = paste0("resources/crohns_case_study/deg/", names(gene_list)[[i]], ".txt"), sep = "\t", col.names = F, quote = F, row.names = F)
}
############################################################################################################################################################
# save annotated MR results - annotate with in DEG (any cell type), DEG_maj_ct, MR_DEG_con, MR_DEG_con_ct_het, 
# mr_annotated <- mr %>% 
#   left_join(deg %>% select(Gene, deg_sig, gut_cell_type, Location, Contrast, `Discrete DE coefficients`, `Discrete DE coefficients p value`, `Discrete FDR`) %>% distinct(), by = join_by(Gene)) %>% 
#   left_join(mr_deg %>% select(Gene, cell_type, deg_maj_ct, mr_deg_sign_concord, mr_deg_sign_ct_het) %>% distinct(), by = join_by(Gene, cell_type)) %>% 
#   mutate(across(where(is.logical), ~replace_na(.x, FALSE)))
# 
# saveRDS(mr_annotated, "resources/crohns_case_study/postprocess/tenk_crohns_sig_deg_appended.RDS")

# prepare supplementary tables 
# annotations <- mr_deg %>%
#   select(Gene, probeID, phenotype, cell_type, major_cell_type, p_SMR_multi, p_HEIDI, gene_type, magma_gene, eqtlgen_mr, canon_gene, drug_gene, drug_pathway_gene, ct_spec, sig_opp_dir_mr_gene, deg_matched_major_ct) %>%
#   distinct()
# dim(annotations)
# 
# annotations_deg <- mr_deg %>%
#   filter(deg_matched_major_ct == TRUE) %>%
#   select(Gene, probeID, cell_type, gut_cell_type, major_cell_type, deg_matched_major_ct, cell_type, concordant_MR_DEG_direction, concordant_MR_DEG_heterogenous_directions, Location, Contrast, `Discrete DE coefficients`, `Discrete DE coefficients p value`, `Discrete FDR`) %>%
#   rename(crohns_dataset_scRNAseq_cell_id = cell_type) %>%
#   distinct()
# 
# dim(annotations_deg)
# unique(annotations_deg$Gene)
# 
# write.csv(annotations, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results.csv", row.names = F)
# saveRDS(annotations, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results.RDS")
# 
# write.csv(annotations_deg, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results_DEG.csv", row.names = F)
# saveRDS(annotations_deg, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results_DEG.RDS")

# # A tibble: 3 × 3
# # Groups:   all_true [2]
# all_true all_false     n
# <lgl>    <lgl>     <int>
#   1 FALSE    FALSE         7
# 2 FALSE    TRUE         47
# 3 TRUE     FALSE        62

# 
# mr_deg_genes <- mr_deg %>% filter(deg_matched_major_ct == T) %>% pull(probeID) %>% unique() 
# concordant_genes <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == T) %>% pull(probeID) %>% unique() 
# concordant_genes_only_true <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == T & concordant_MR_DEG_heterogenous_directions == 1) %>% pull(probeID) %>% unique() 
# discordant_genes <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == F) %>% pull(probeID) %>% unique()
# discordant_genes_only_false <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == F & concordant_MR_DEG_heterogenous_directions == 1) %>% pull(probeID) %>% unique()

# note only true and only false refer to whether its only concordant or only discordant, because some genes are both.
#save(mr_deg_genes, concordant_genes, concordant_genes_only_true, discordant_genes, discordant_genes_only_false, file = "resources/crohns_case_study/crohns_relevant_gene_lists/mr_deg_gene_lists.RData")
