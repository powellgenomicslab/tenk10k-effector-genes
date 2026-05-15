# qsub -I -q normal -P fy54 -l ncpus=4,storage=gdata/fy54+gdata/ei56,mem=50GB -l jobfs=100GB
tenk_dir <- paste0("/g/data/fy54/analysis/tenk10k-causal")
setwd(tenk_dir)
crohns_dir <- "resources/crohns_case_study/postprocess"

# source preprocess scripts 
source(paste0(tenk_dir, "/scripts/preprocess.R"))

# save results from preprocess script
saveRDS(df_msmr_tenk10k, paste0(crohns_dir, "/tenk_alltraits_all.RDS"))
saveRDS(df_msmr, paste0(crohns_dir, "/tenk_alltraits_sig.RDS"))

# filter to get crohns only all and sig results 
tenk_crohns_all <- filter(df_msmr_tenk10k, phenotype == "crohns")
tenk_crohns_sig <- filter(df_msmr, phenotype == "crohns")

saveRDS(tenk_crohns_all, paste0(crohns_dir,"/tenk_crohns_all.RDS"))
saveRDS(tenk_crohns_sig, paste0(crohns_dir,"/tenk_crohns_sig.RDS"))

# get the Gene names
get_annotation <- function(df, gtf_df) {
  df <- df %>%
    left_join(gtf_df %>% select(hgnc_symbol, ensembl_gene_id), by = c("probeID" = "ensembl_gene_id"))
  return(df)
}

# change Gene column from ens ids to gene names
df_msmr_eqtlgen <- get_annotation(df_msmr_eqtlgen, df_gene_annot)
df_msmr_eqtlgen <- df_msmr_eqtlgen %>%
  select(-Gene) %>%
  rename(Gene = hgnc_symbol)

# save eqtlgen
saveRDS(df_msmr_eqtlgen, paste0(crohns_dir, "/eqtlgen_alltraits_all.RDS"))

# filter to sig eqtlgen all traits and sig and all for crohns only 
eqtlgen_alltraits_sig <- filter(df_msmr_eqtlgen, sig == TRUE)
eqtlgen_crohns_all <- filter(df_msmr_eqtlgen, phenotype == "crohns")
eqtlgen_crohns_sig <- filter(eqtlgen_alltraits_sig, phenotype == "crohns")

saveRDS(eqtlgen_alltraits_sig, paste0(crohns_dir,"/eqtlgen_alltraits_sig.RDS"))
saveRDS(eqtlgen_crohns_all, paste0(crohns_dir,"/eqtlgen_crohns_all.RDS"))
saveRDS(eqtlgen_crohns_sig, paste0(crohns_dir,"/eqtlgen_crohns_sig.RDS"))


