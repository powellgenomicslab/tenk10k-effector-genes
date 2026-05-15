library(tidyverse)

#crohns_dir = "/g/data/fy54/analysis/tenk10k-causal/resources/crohns_case_study"
crohns_dir = "resources/crohns_case_study"

mr <- readRDS(paste0(crohns_dir, "/postprocess/tenk_crohns_sig.RDS")) 
mr_all <- readRDS(paste0(crohns_dir, "/postprocess/tenk_crohns_all.RDS")) 

# how many genes with opposite directions of effect in different cell types, and how many of those protein coding 

# create column to annotate with opposite directions of effect. 
mr <- mr %>% 
  group_by(Gene) %>% 
  mutate(across(c(b_GWAS, b_SMR, b_eQTL), ~var(sign(.x)) == 0, .names = "concordant.{.col}")) 
  # %>% 
  # filter(gene_type == "protein_coding")

mr_all <- mr_all %>% 
  group_by(Gene) %>% 
  mutate(across(c(b_GWAS, b_SMR, b_eQTL), ~var(sign(.x)) == 0, .names = "concordant.{.col}")) 
  # %>% 
  # filter(gene_type == "protein_coding")

# How many unique genes have concordant bSMR effect directions
opp_dir_df <- mr %>% group_by(concordant.b_SMR) %>% summarize(genes = n_distinct(Gene))
# NA column corresponds to cell type specific entries where there is no heterogeneity

# # A tibble: 3 x 2
#   concordant.b_SMR genes
#   <lgl>            <int>
# 1 FALSE               61
# 2 TRUE               243
# 3 NA                 206

saveRDS(mr, paste0(crohns_dir, "/postprocess/tenk_crohns_sig.RDS"))
saveRDS(mr_all, paste0(crohns_dir, "/postprocess/tenk_crohns_all.RDS"))

# which genes have opposite directions of effect driven by eqtl differences not gwas differences 
# either gwas concordant or eqtl discordant 

mr_eqtl_discordant <- mr %>% filter(concordant.b_SMR == F, concordant.b_eQTL == F) 
mr_gwas_concordant <- mr %>% filter(concordant.b_SMR == F, concordant.b_GWAS) 
mr_smr_disconcordant <- mr %>% filter(concordant.b_SMR == F)

opp_dir_genes <- unique(mr_smr_disconcordant$probeID)
opp_dir_genes_concordant_bGWAS <- unique(mr_gwas_concordant$probeID) 
opp_dir_genes_discordant_beqtl <- unique(mr_eqtl_discordant$probeID) 

length(opp_dir_genes) # 61
length(opp_dir_genes_concordant_bGWAS) # 20
length(opp_dir_genes_discordant_beqtl ) # 43

#################################################################################################

# write.csv(mr, "causal_inference_manuscript/supplementary/tenk_crohns_sig_annotate_directions_bGWAS_beQTL_bSMR.csv", row.names = F)
# write.csv(mr_gwas_concordant, "causal_inference_manuscript/supplementary_reports/tenk_crohns_sig_annotate_directions_concordant_bGWAS_disconcordant_bSMR.csv", row.names = F)
#saveRDS(opp_dir_genes_concordant_bGWAS, "causal_inference_manuscript/data/crohns_relevant_gene_lists/opp_dir_genes_concordant_bGWAS_sig.RDS")

# Save the significant opposite direction genes, for simply opposite bMR but also opposite bMR + opposite beQTL

write.csv(mr_smr_disconcordant, paste0(crohns_dir, "/supplementary/tenk_crohns_sig_annotate_directions_discordant_bSMR.csv"), row.names = F)
saveRDS(opp_dir_genes, paste0(crohns_dir, "/crohns_relevant_gene_lists/sig_opp_dir_genes.RDS"))

write.csv(mr_eqtl_discordant, paste0(crohns_dir, "/supplementary/tenk_crohns_sig_annotate_directions_discordant_bSMR_discordant_eQTL.csv"), row.names = F)
saveRDS(opp_dir_genes_discordant_beqtl, paste0(crohns_dir, "/crohns_relevant_gene_lists/sig_opp_dir_genes_discordant_beqtl.RDS"))
#################################################################################################

