library(tidyverse)

crohns_dir = "/g/data/fy54/analysis/tenk10k-causal/resources/crohns_case_study"

# the table X came from manually extracting names from two papers Tyebally et al and __ et al 
X <- read.table(paste0(crohns_dir, "/crohns_relevant_gene_lists/sources/drug_target_genes.txt"))
drug_target_genes <- X$V1
saveRDS(drug_target_genes, paste0(crohns_dir, "/crohns_relevant_gene_lists/drug_target_genes.RDS"))

# these lists are copied from sources, check sources folder
tnf <- c("MAP2K7,IKBKB,CHUK,MAP3K7,SMPD2,BAG4,MAP4K4,TNF,TXN,SMPD1,TNFRSF1A,NFKB1,TNFRSF1B,TNFAIP3,PRKCI,STAT1,MAP2K3,RACK1,ADAM17,CAV1,RELA,PRKCZ,MAP4K2,TRAF2,TRAF1,FADD,MAP3K1,BIRC3,BIRC2,SQSTM1,RIPK1,CASP8,TRADD,TAB1,NRK,MAP4K3,MADD,RFFL,NSMAF,MAP3K5,MAP3K3,CYLD,TAB2,TNIK,MAP4K5,IKBKG")
tnf <- str_split_1(tnf, ",")

il12_il23 <- c("SOCS1,RIPK2,JAK2,GADD45B,IL18RAP,GADD45G,EOMES,FOS,IFNG,IL1B,IL2RA,CD4,CD8A,HLA-DRA,,CD3D,HLA-A,IL4,LCK,CD3E,CD3G,GZMB,CCL3,CD8B,GZMA,CCL4,IL1R1,IL2RB,ATF2,PPP3CB,NFKB1,CD247,IL12A,IL12B,TYK2,IL2RG,NOS2,STAT3,STAT1,STAT6,STAT5A,MTOR,IL12RB1,MAP2K3,FASLG,RAB7A,CCR5,MAP2K6,IL2,B2M,PPP3R1,NFKB2,RELB,RELA,PPP3CA,IL18R1,IL18,STAT4,HLX,MAPK14,IL12RB2,SPHK2,TBX21,SOCS3,JAK2,ALOX12B,IL18RAP,TNF,IFNG,IL1B,CD4,MPO,IL6,CD3E,CXCL1,CCL2,NFKB1,NFKBIA,ITGA3,PIK3R1,IL12B,TYK2,NOS2,STAT3,STAT1,STAT5A,PIK3CA,IL12RB1,IL2,RELA,CXCL9,IL24,IL18R1,IL18,STAT4,IL17A,IL23R,IL17F,IL23A,IL19")
il12_il23 <- str_split_1(il12_il23, ",")
il12_il23 <- unique(il12_il23)

integrin <- c("ITGA10,ITGB3,ITGB2,ITGB1,ITGAV,ITGA2B,ITGA5,ITGAM,ITGA4,ITGB4,ITGA2,ITGB5,ITGB6,ITGAL,ITGAX,ITGA6,ITGA3,ITGB7,ITGB8,ITGAE,ITGA8,ITGA1,ITGAD,ITGA7,ITGA9,ITGA11,LAMA5,ITGA10,F13A1,PLAU,COL1A1,COL2A1,COL3A1,COL4A1,FGA,FGB,FGG,FN1,VTN,ITGB1,COL5A2,ITGAV,LAMB1,THBS1,COL1A2,CD14,ITGA5,SPP1,LAMC1,COL11A1,COL6A1,COL6A2,COL6A3,ITGA4,COL11A2,NID1,VEGFA,ITGA2,VCAM1,COL5A1,MDK,TGM2,ITGA6,LAMA2,TNC,LAMA1,ITGA3,COL4A5,THBS2,FBN1,COL18A1,COL4A4,ITGA8,LAMB2,ITGA1,JAM2,CD81,COL4A3,COL7A1,PLAUR,ITGA7,LAMB3,LAMC2,ITGA9,COL4A6,TGFBI,LAMA4,LAMA3,CSPG4,NPNT,IGSF8,ITGA11,CCN1,EDIL3,PLAU,FN1,VTN,ITGAV,ITGA4,ITGB5,ITGB6,SDC1,VCAM1,ITGB7,ITGB8,FBN1,TGFBR1,PLAUR,MADCAM1")
integrin <- str_split_1(integrin, ",")
integrin <- unique(integrin)

all_drug_pathway_genes <- unique(c(tnf, il12_il23, integrin, drug_target_genes))

saveRDS(all_drug_pathway_genes, paste0(crohns_dir, "/crohns_relevant_gene_lists/drug_pathway_genes.RDS"))
# saveRDS(il12_il23, "causal_inference_manuscript/data/crohns_relevant_gene_lists/il12_il23_genes.RDS")
# saveRDS(integrin, "causal_inference_manuscript/data/crohns_relevant_gene_lists/integrin_genes.RDS")
# saveRDS(tnf, "causal_inference_manuscript/data/crohns_relevant_gene_lists/tnf_genes.RDS")

save.image(paste0(crohns_dir, "/crohns_relevant_gene_lists/prepare_drug_pathway_genes.RData"))

