
library(tidyverse)
#crohns_dir = "/g/data/fy54/analysis/tenk10k-causal/resources/crohns_case_study"
crohns_dir = "resources/crohns_case_study"
canon_genes <- read.csv(paste0(crohns_dir, "/crohns_relevant_gene_lists/sources/Systematic Review of Crohns Disease Genes.csv"))

canon_genes_categories <- unique(canon_genes$Category)

# [1] "Experimental evidence of variant"                 "Other evidences of genetic alterations"           "Treatment Response"                              
# [4] "Annotation error"                                 "Biologically related but no evidence of mutation" "Genetic evidence in a related disease"           
# [7] "GWAS evidence within gene"                        "Unrelated"                                        "Non-Human"                                       
# [10] "Genetic evidence in related complications"        "Negative evidence"  

# select a few only "Experimental evidence of variant","Other evidences of genetic alterations", "Treatment Response" , "Biologically related but no evidence of mutation" 
# "GWAS evidence within gene"     
canon_genes_categories <- canon_genes_categories[c(1,2,3,5,7)]

canon_genes <- canon_genes %>% 
  filter(Category %in% canon_genes_categories) %>% 
  filter(Lvl >=8, Abstracts >= 10) %>%  # Select Document score above 8 and number of Abstracts above 10
  arrange(desc(Abstracts), Lvl) %>% 
  select(-X.)

write.csv(canon_genes, paste0(crohns_dir, "/supplementary/canonical_gene_info_table_etal_2023.csv"), row.names = F)

canon_genes_genes <- canon_genes %>% pull(Gene)

saveRDS(canon_genes_genes, paste0(crohns_dir, "/crohns_relevant_gene_lists/canonical_genes.RDS"))
