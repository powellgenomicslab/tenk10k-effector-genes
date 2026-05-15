# Merge DEG and MR by major cell type and save .csv and gene list 

library(tidyverse) 

deg <- readRDS("resources/crohns_case_study/deg/crohns_deg_pre-processed.RDS") %>% 
  select(Gene, Contrast, scRNAseq_cellid, everything())

mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% 
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR)) 

cell_map <- read.delim("resources/metadata/cell_map.tsv") %>% as.data.frame()

# check
any(is.na(mr$p_transform))

# filter deg genes
deg <- deg %>% filter(Gene %in% mr$Gene)

############################################################################################################################################################

# Merged the MR and DEG dataframes on gene and major cell type
mr_deg <- left_join(mr, deg, by = c("Gene", "major_cell_type")) %>% 
  mutate(deg_matched_major_ct = !is.na(Contrast)) %>% 
  mutate(concordant_MR_DEG_direction = sign(b_SMR) == sign(`Discrete DE coefficients`)) %>% 
  group_by(Gene) %>% 
  mutate(concordant_MR_DEG_heterogenous_directions = n_distinct(concordant_MR_DEG_direction)) %>% 
  ungroup()

# save to deg
saveRDS(mr_deg, "resources/crohns_case_study/deg/MR_and_DEG_matched_major_cell_type_combined_results_colon_df.RDS")

############################################################################################################################################################
# prepare supplementary tables 

annotations <- mr_deg %>%
  select(Gene, probeID, phenotype, cell_type, major_cell_type, p_SMR_multi, p_HEIDI, gene_type, magma_gene, eqtlgen_mr, canon_gene, drug_gene, drug_pathway_gene, ct_spec, sig_opp_dir_mr_gene, deg_matched_major_ct) %>% 
  distinct()
dim(annotations)

annotations_deg <- mr_deg %>%
  filter(deg_matched_major_ct == TRUE) %>% 
  select(Gene, probeID, cell_type, major_cell_type, deg_matched_major_ct, scRNAseq_cellid, concordant_MR_DEG_direction, concordant_MR_DEG_heterogenous_directions, Location, Contrast, `Discrete DE coefficients`, `Discrete DE coefficients p value`, `Discrete FDR`) %>% 
  rename(crohns_dataset_scRNAseq_cell_id = scRNAseq_cellid) %>%  
  distinct()

dim(annotations_deg)
unique(annotations_deg$Gene)

write.csv(annotations, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results.csv", row.names = F)
saveRDS(annotations, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results.RDS")

write.csv(annotations_deg, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results_DEG.csv", row.names = F)
saveRDS(annotations_deg, "resources/crohns_case_study/supplementary/Crohns_Annotated_Results_DEG.RDS")


# how many genes have concordant directions of effect (true or false), how many genes have both concordant and discordant effects depending on cell types?
number_of_genes_grouped <- mr_deg %>% 
  filter(deg_matched_major_ct == T) %>% 
  group_by(Gene) %>% 
  summarize(all_true = all(concordant_MR_DEG_direction), all_false = all(!concordant_MR_DEG_direction)) %>% 
  group_by(all_true, all_false) %>% 
  tally() 
# 
# # A tibble: 3 × 3
# # Groups:   all_true [2]
# all_true all_false     n
# <lgl>    <lgl>     <int>
#   1 FALSE    FALSE         7
# 2 FALSE    TRUE         47
# 3 TRUE     FALSE        62


mr_deg_genes <- mr_deg %>% filter(deg_matched_major_ct == T) %>% pull(probeID) %>% unique() 
concordant_genes <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == T) %>% pull(probeID) %>% unique() 
concordant_genes_only_true <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == T & concordant_MR_DEG_heterogenous_directions == 1) %>% pull(probeID) %>% unique() 
discordant_genes <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == F) %>% pull(probeID) %>% unique()
discordant_genes_only_false <- mr_deg %>% filter(deg_matched_major_ct == T & concordant_MR_DEG_direction == F & concordant_MR_DEG_heterogenous_directions == 1) %>% pull(probeID) %>% unique()

# note only true and only false refer to whether its only concordant or only discordant, because some genes are both.
save(mr_deg_genes, concordant_genes, concordant_genes_only_true, discordant_genes, discordant_genes_only_false, file = "resources/crohns_case_study/crohns_relevant_gene_lists/mr_deg_gene_lists.RData")
