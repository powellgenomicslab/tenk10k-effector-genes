# Date: 07/04/2026 
# Purpose: Prepare external source annotation dataframe 
# 1: Canonical genes from Garza-Hernandez et al 2022 
# 2: Any evidence of association with Crohn's disease on OpenTargetsPlatform 
# 3: Nearest gene to 320 loci from Liu et al., 2023 
# 4: Drug target list from literature 
# 5: Drug pathway list from GSEA pathways 

library(here)
library(tidyverse)

# Gene annotation table 
df_gene_annot <- data.table::fread("resources/misc/gencode.v44.gene_type.tsv")

############# For every source, generate a dataframe with the column name "Gene"############

############ 1: Canonical genes from Garza-Hernandez et al 2022 PMID: 35418025
# Read in canonical genes 
df1 <- read.csv(here("resources/crohns_case_study/crohns_relevant_gene_lists/sources/Systematic Review of Crohns Disease Genes.csv"), row.names = 1)
df1_categories <- unique(df1$Category)
# print(df1_categories)
#  [1] "Experimental evidence of variant"                
#  [2] "Other evidences of genetic alterations"          
#  [3] "Treatment Response"                              
#  [4] "Annotation error"                                
#  [5] "Biologically related but no evidence of mutation"
#  [6] "Genetic evidence in a related disease"           
#  [7] "GWAS evidence within gene"                       
#  [8] "Unrelated"                                       
#  [9] "Non-Human"                                       
# [10] "Genetic evidence in related complications"       
# [11] "Negative evidence"   
df1_categories <- df1_categories[c(1,2,3,5,7)] # select:  "Experimental evidence of variant", "Other evidences of genetic alterations", "Treatment Response" , "Biologically related but no evidence of mutation", "GWAS evidence within gene"     

# Filter the Systematic review results to high confidence ones to call as "canonical gene list"
df1 <- df1 %>% 
    as_tibble() %>%
    filter(Category %in% df1_categories) %>% 
    # Select Document score above 8 and number of Abstracts above 10
    filter(Lvl >=8, Abstracts >= 10) %>%  
    select(Gene) %>%
    distinct()
############ 2: Any evidence of association with Crohn's disease on OpenTargetsPlatform 
# Note: Only Clinical and Genetic Associations exported from: https://platform.opentargets.org/disease/EFO_0000384/associations
df2 <- read_tsv(here("resources/crohns_case_study/crohns_relevant_gene_lists/sources/OT-EFO_0000384-associated-targets-07_04_2026-v26_03.tsv"), show_col_types = FALSE) 
# Global score 
df2 <- df2 %>% dplyr::rename(Gene = symbol) %>% 
  dplyr::filter(globalScore > quantile(df2$globalScore, probs = 0.95)) %>% 
  dplyr::select(Gene) %>% 
  dplyr::distinct()

############ 3: Nearest gene to 320 loci from Liu et al., 2023 
df3 <- readxl::read_xlsx(here("resources/crohns_case_study/crohns_relevant_gene_lists/sources/41588_2023_1384_MOESM4_ESM_Liuetal2023_SupTable.xlsx"), sheet = "ST8", skip = 1) %>% select(Gene) %>% distinct()

############# 4: Drug target list from literature 
# Names of drug targets from the literature
# PMID: 36587559, Table 1. 
# PMID: 37983917, Table 1. 
df4 <- readLines(here("resources/crohns_case_study/crohns_relevant_gene_lists/sources/drug_target_genes.txt")) %>% as_tibble_col(column_name = "Gene") %>% distinct()

############# 5: Drug pathway list from GSEA pathway lists
df5a <- c("MAP2K7,IKBKB,CHUK,MAP3K7,SMPD2,BAG4,MAP4K4,TNF,TXN,SMPD1,TNFRSF1A,NFKB1,TNFRSF1B,TNFAIP3,PRKCI,STAT1,MAP2K3,RACK1,ADAM17,CAV1,RELA,PRKCZ,MAP4K2,TRAF2,TRAF1,FADD,MAP3K1,BIRC3,BIRC2,SQSTM1,RIPK1,CASP8,TRADD,TAB1,NRK,MAP4K3,MADD,RFFL,NSMAF,MAP3K5,MAP3K3,CYLD,TAB2,TNIK,MAP4K5,IKBKG")
df5a <- str_split_1(df5a, ",") %>% as_tibble_col(column_name = "Gene") %>% distinct()

df5b <- c("SOCS1,RIPK2,JAK2,GADD45B,IL18RAP,GADD45G,EOMES,FOS,IFNG,IL1B,IL2RA,CD4,CD8A,HLA-DRA,,CD3D,HLA-A,IL4,LCK,CD3E,CD3G,GZMB,CCL3,CD8B,GZMA,CCL4,IL1R1,IL2RB,ATF2,PPP3CB,NFKB1,CD247,IL12A,IL12B,TYK2,IL2RG,NOS2,STAT3,STAT1,STAT6,STAT5A,MTOR,IL12RB1,MAP2K3,FASLG,RAB7A,CCR5,MAP2K6,IL2,B2M,PPP3R1,NFKB2,RELB,RELA,PPP3CA,IL18R1,IL18,STAT4,HLX,MAPK14,IL12RB2,SPHK2,TBX21,SOCS3,JAK2,ALOX12B,IL18RAP,TNF,IFNG,IL1B,CD4,MPO,IL6,CD3E,CXCL1,CCL2,NFKB1,NFKBIA,ITGA3,PIK3R1,IL12B,TYK2,NOS2,STAT3,STAT1,STAT5A,PIK3CA,IL12RB1,IL2,RELA,CXCL9,IL24,IL18R1,IL18,STAT4,IL17A,IL23R,IL17F,IL23A,IL19")
df5b <- str_split_1(df5b, ",") %>% as_tibble_col(column_name = "Gene") %>% distinct()

df5c <- c("ITGA10,ITGB3,ITGB2,ITGB1,ITGAV,ITGA2B,ITGA5,ITGAM,ITGA4,ITGB4,ITGA2,ITGB5,ITGB6,ITGAL,ITGAX,ITGA6,ITGA3,ITGB7,ITGB8,ITGAE,ITGA8,ITGA1,ITGAD,ITGA7,ITGA9,ITGA11,LAMA5,ITGA10,F13A1,PLAU,COL1A1,COL2A1,COL3A1,COL4A1,FGA,FGB,FGG,FN1,VTN,ITGB1,COL5A2,ITGAV,LAMB1,THBS1,COL1A2,CD14,ITGA5,SPP1,LAMC1,COL11A1,COL6A1,COL6A2,COL6A3,ITGA4,COL11A2,NID1,VEGFA,ITGA2,VCAM1,COL5A1,MDK,TGM2,ITGA6,LAMA2,TNC,LAMA1,ITGA3,COL4A5,THBS2,FBN1,COL18A1,COL4A4,ITGA8,LAMB2,ITGA1,JAM2,CD81,COL4A3,COL7A1,PLAUR,ITGA7,LAMB3,LAMC2,ITGA9,COL4A6,TGFBI,LAMA4,LAMA3,CSPG4,NPNT,IGSF8,ITGA11,CCN1,EDIL3,PLAU,FN1,VTN,ITGAV,ITGA4,ITGB5,ITGB6,SDC1,VCAM1,ITGB7,ITGB8,FBN1,TGFBR1,PLAUR,MADCAM1")
df5c <- str_split_1(df5c, ",") %>% as_tibble_col(column_name = "Gene") %>% distinct()

df6 <- dplyr::bind_rows(df1, df2, df3, df4, df5a, df5b, df5c) %>% distinct()
############# Add more gene lists above ############

############# Create a list of sources #############

gene_lists <- list(
  sys_review = df1, 
  open_targets_cd = df2, 
  liuetal = df3,
  drug_target = df4, 
  drug_pathway_tnf = df5a,
  drug_pathway_il23 = df5b, 
  drug_pathway_integrin = df5c,
  cd_known = df6
)
# Create a dataframe of sources and genes and a membership status column (placeholder for wide format next step)
gene_lists_df <- dplyr::bind_rows(gene_lists, .id = "source") %>% mutate(member = TRUE) %>% distinct()

# Pivot wider to create annotation table 
gene_lists_df <- gene_lists_df %>% 
    pivot_wider(names_from = source, values_from = member, values_fill = FALSE) %>%
    left_join(df_gene_annot %>% select(hgnc_symbol, ensembl_gene_id), by = c("Gene" = "hgnc_symbol")) %>%
    select(Gene, ensembl_gene_id, everything()) %>%
    drop_na()

length(unique(gene_lists_df$Gene))
head(gene_lists_df)

saveRDS(gene_lists_df, here("resources/crohns_case_study/crohns_relevant_gene_lists/external_source_gene_annot_df.RDS"))