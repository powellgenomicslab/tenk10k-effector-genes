# get sig canonical and drug pathway genes. this is done post annotation 

library(tidyverse)
crohns_dir = "/g/data/fy54/analysis/tenk10k-causal/resources/crohns_case_study"

res <- readRDS(paste0(crohns_dir, "/postprocess/tenk_crohns_sig.RDS"))

canon_genes <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/canonical_genes.RDS"))

sig_canon_genes <- res %>% 
  filter(canon_gene == T) %>% 
  filter(sig == T) %>% 
  pull(probeID)

sig_drug_pathway_genes <- res %>% 
  filter(drug_gene == T | drug_pathway_gene == T) %>% 
  filter(sig == T) %>% 
  pull(probeID)

saveRDS(sig_canon_genes, paste0(crohns_dir, "/crohns_relevant_gene_lists/sig_canon_genes.RDS"))
saveRDS(sig_drug_pathway_genes, paste0(crohns_dir, "/crohns_relevant_gene_lists/sig_drug_pathway_or_target_genes.RDS"))
