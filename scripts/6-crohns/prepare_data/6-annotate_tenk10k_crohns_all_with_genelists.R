library(tidyverse)
library(data.table)

crohns_dir = "/g/data/fy54/analysis/tenk10k-causal/resources/crohns_case_study"

# load all the relevant lists 

canon <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/canonical_genes.RDS"))
drug <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/drug_target_genes.RDS"))
drug_pathway <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/drug_pathway_genes.RDS"))
eqtlgen <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/eqtlgen_crohns_sig_genes.RDS"))
ct_spec_genes <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/tenk_crohns_sig_cell_type_specific_genes.RDS"))
sig_opp_dir_gene <- readRDS(paste0(crohns_dir, "/crohns_relevant_gene_lists/sig_opp_dir_genes.RDS"))

# add respective annotation columns 
# the eqtlgen_sig is redundant now, as the new preprocess script has this, but leaving it for harmony with other script

annotate_df <- function(cd){
  cd$canon_gene <- ifelse(cd$Gene %in% canon , T, F)
  cd$drug_gene <- ifelse(cd$Gene %in% drug , T, F)
  cd$drug_pathway_gene <- ifelse(cd$Gene %in% drug_pathway, T, F)
  #cd$eqtlgen_sig <- ifelse(cd$probeID %in% eqtlgen$probeID , T, F)
  cd$ct_spec <- ifelse(cd$probeID %in% ct_spec_genes, T, F)
  cd$sig_opp_dir_mr_gene <- ifelse(cd$probeID %in% sig_opp_dir_gene, T, F)
  return(cd)
}

#list.files(paste0(crohns_dir, "/crohns_post_process/")[c(3,4)]
all <- readRDS(paste0(crohns_dir, "/postprocess/tenk_crohns_all.RDS"))
sig <- readRDS(paste0(crohns_dir, "/postprocess/tenk_crohns_sig.RDS"))

all <- annotate_df(all)
sig <- annotate_df(sig)

# forgot to hash out of the function earlier
# sig <- sig |> select(-eqtlgen_sig)
# all <- all |> select(-eqtlgen_sig)
saveRDS(all, paste0(crohns_dir, "/postprocess/tenk_crohns_all.RDS"))
saveRDS(sig, paste0(crohns_dir, "/postprocess/tenk_crohns_sig.RDS"))
