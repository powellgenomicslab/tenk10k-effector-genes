annot <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/external_source_gene_annot_df.RDS") %>% 
  rename("Systematic Review" = sys_review,
         "95th Pctl Open Targets" = open_targets_cd,
         "IBD metaGWAS Nearest Gene" = liuetal,
         "Drug Target" = drug_target,
         "TNF Pathway" = drug_pathway_tnf,
         "IL23 Pathway" = drug_pathway_il23, 
         "Integrin Pathway" = drug_pathway_integrin,
         "Any Annotation" = cd_known,
         "EnsemblID" = ensembl_gene_id)

