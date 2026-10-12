# Prepare supplementary table - all CD genes in other immune diseases.

library(tidyverse)

mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% 
  mutate(across(where(is.logical), function(x) (replace_na(x, FALSE))))

#mr <- mr %>% filter(ct_spec == TRUE & eqtlgen_sig == FALSE)
#mr <- mr %>% filter(eqtlgen_sig == FALSE)
goi <- mr %>% pull(Gene) %>% unique() 

mr_pheno_ct <- readRDS("resources/crohns_case_study/postprocess/tenk_alltraits_sig.RDS") %>% 
  filter(supercategory == "disease", pheno_cat == "Immune") %>% 
  filter(Gene %in% goi) %>% 
  select(pheno_label, Gene, cell_type, probeID) %>% 
  group_by(Gene, probeID) %>% 
  summarise(Phenotypes = paste(unique(pheno_label), collapse = "; "), 
            `Cell Types` = paste(unique(cell_type), collapse = "; "))

# Get the annotations in the same way 
annot_cols <- c("eqtlgen_sig", "sys_review", "open_targets_cd", "liuetal", "drug_target", "drug_pathway_tnf", 
                "drug_pathway_il23", "drug_pathway_integrin")

rename_cols <- c("Gene", "probeID", "Bulk eQTLGen MR", "Systematic Review", "95th Pctl Open Targets", "IBD metaGWAS Nearest Gene", "Drug Target", "TNF Pathway", 
                 "IL23 Pathway", "Integrin Pathway")

mr_annot <- mr %>% 
  filter(Gene %in% goi) %>% 
  dplyr::select(Gene, probeID,all_of(annot_cols)) %>% 
  rename_with(~ rename_cols) %>% 
  pivot_longer(cols = c(-Gene,-probeID), names_to = "annotation_cols", values_to = "Present") %>%
  mutate(annotation = ifelse(Present, annotation_cols, NA)) %>%
  group_by(Gene, probeID) %>% 
  summarize(`Crohn's Disease Gene Annotation` = paste(unique(na.omit(annotation)), collapse = "; "))

crohns_immune <- left_join(mr_pheno_ct, mr_annot, by = c("Gene", "probeID"))

#saveRDS(res, "resources/crohns_case_study/supplementary/cd_genes_other_diseases_and_annotations_sup_table.RDS")
