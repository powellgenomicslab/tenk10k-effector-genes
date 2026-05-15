# Crohn's MR association summary statistics, with gene annotations
annot <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/external_source_gene_annot_df.RDS")
annot_names <- colnames(annot)[!colnames(annot) == c("Gene", "ensembl_gene_id")]
mr_deg <- readRDS("resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS") 

crohns_annotated_results <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% 
  mutate(major_cell_type = as.character(major_cell_type) %>% dplyr::replace_when(cell_type == "ILC" ~ "NK", cell_type == "Plasmablast" ~ "Plasma B")) %>%
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR)) %>% 
  left_join(mr_deg %>% select(Gene, major_cell_type, deg_maj_ct), by = join_by(Gene, major_cell_type)) %>% 
  select(Gene, probeID, cell_type, b_SMR, p_transform, p_SMR_multi, p_HEIDI, phet_ivw, psigmay_mrlink2, coloc_pph4, mvcoloc_pph4, discordant_bSMR, magma_gene, eqtlgen_sig, all_of(annot_names), deg_maj_ct) %>% 
  mutate(across(where(is.logical), function(x) replace_na(x, FALSE))) %>% 
  rename("Cell Type" = cell_type, 
         "EnsemblID" = probeID,
         "Single SNP beta MR" = b_SMR,
         "Signed log10 P MR" = p_transform,
         `P Multi-SNP MR` = p_SMR_multi, 
         `P HEIDI` = p_HEIDI,
         `Cochran's Q` = phet_ivw, 
         `P Σy MRLink2` = psigmay_mrlink2, 
         `Single-variant COLOC PP H4` = coloc_pph4, 
         `Multi-variant COLOC PP H4` = mvcoloc_pph4, 
         "MAGMA" = magma_gene,
         "Bulk eQTLGen MR" = eqtlgen_sig, 
         "Systematic Review" = sys_review,
         "95th Pctl Open Targets" = open_targets_cd,
         "IBD metaGWAS Nearest Gene" = liuetal,
         "Drug Target" = drug_target,
         "TNF Pathway" = drug_pathway_tnf,
         "IL23 Pathway" = drug_pathway_il23, 
         "Integrin Pathway" = drug_pathway_integrin,
         "Any Literature Annotation" = cd_known,
         "Discordant MR Effect Between Cell Types" = discordant_bSMR,
         "Differentially Expressed in Matched Major Cell Type (Kong et al., 2023)" = deg_maj_ct) 
