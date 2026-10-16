# Purpose:
# Create a major cell type classification for Kong et al., 2023 DEG results 
# Create an updated major cell type classification for Cuomo et al., 2025 from metadata/cell.tsv
# Inner-join MR and DEG results on major cell type classification 

#########################################################################################
# Interactive job information 
# qsub -I -q normal -P ei56 -l ncpus=2,storage=gdata/fy54+gdata/ei56,mem=8GB -l jobfs=100GB
# conda activate tidyverse
# R
#########################################################################################

library(tidyverse)
library(here) 

crohns_dir <- here("resources/crohns_case_study")
# Kong et al. 2023 DEG table, immune cells and discrete model only, with Ensembl IDs
# (extracted from the published supplementary table by
# workflow/rules/snakescripts/crohns_external/extract_third_party.R)
deg <- arrow::read_parquet(paste0(crohns_dir, "/deg/kong2023_immune_deg.parquet"))

print("Kong et al DEG colnames: ")
colnames(deg)

deg <- deg %>% rename(gut_cell_type = `Cell subset`)
deg$gut_cell_type <- gsub("Immune.", "", as.character(deg$gut_cell_type))

# Prepare a -log10 transformed pvalue that also shows the direction of effect by multiplying with direction of beta (-log10 Discrete pval * direction of beta)
deg$disc_p_transform <- -log10(deg$`Discrete DE coefficients p value`)*sign(deg$`Discrete DE coefficients`) 

# Convert location (either colon or terminal ileum) and contrast (either inflamed vs healthy or healthy vs inflamed) to factors with levels 
# deg <- deg %>% mutate(Location = ifelse(Location == "CO", "Colon", "Terminal Ileum"), Contrast = ifelse(Contrast == "Infl vs. Heal", "Inflamed Tissue vs Healthy", "Non-Inflamed Tissue vs Healthy"))
deg$Location <- factor(deg$Location, levels = c("CO", "TI"))
deg$Contrast <- factor(deg$Contrast, levels = c("Infl vs. Heal", "NonI vs. Heal"))


print("Explore scRNAseq ID cell types: ")
deg_celltypes <- unique(deg$gut_cell_type) |> sort()
print(deg_celltypes)

# create major cell type annotation 
maj_annot <- deg %>% 
  select(gut_cell_type) %>% 
  distinct() %>% 
  mutate(major_cell_type = case_when(
    gut_cell_type %in% c("T cells CD4+ FOSB+", "T cells CD4+ IL17A+","T cells Naive CD4+","Tregs") ~ "CD4 T", 
    gut_cell_type %in% c("T cells CD8+","T cells CD8+ KLRG1+") ~ "CD8 T", 
    gut_cell_type %in% c("T cells OGT+") ~ "Unconventional T",
    gut_cell_type %in% c("Plasma cells" ) ~ "Plasma B",
    gut_cell_type %in% c("NK cells KLRF1+ CD3G-", "ILCs") ~ "NK", 
    gut_cell_type %in% c("B cells") ~ "B",
    gut_cell_type %in% c("Monocytes S100A8+ S100A9+") ~ "Monocyte",
    gut_cell_type %in% c("DC1","DC2 CD1D+","DC2 CD1D-","Mature DCs") ~ "Dendritic",
    gut_cell_type %in% c("Mast cells") ~ "Mast",
    gut_cell_type %in% c("Macrophages", "Macrophages CCL3+ CCL4+","Macrophages CXCL9+ CXCL10+","Macrophages LYVE1+","Macrophages Metallothionein","Macrophages PLA2G2D+") ~ "Macrophage",
    gut_cell_type %in% c("NK-like cells ID3+ ENTPD1+", "IELs ID3+ ENTPD1+") ~ "Intra-epithelial Lymphocytes",
    gut_cell_type %in% c("Cycling cells") ~ "Cycling",
    gut_cell_type %in% c("B cells AICDA+ LRMP+") ~ "Germinal Centre B")) %>%
  mutate(tissue_specific = if_else(!major_cell_type %in% c("CD4 T", "CD8 T", "Unconventional T", "B", "Plasma B", "NK", "Dendritic", "Monocyte"), TRUE, FALSE))

# Create a major cell type annotation column in the DEG results to match with MR results later 
deg <- deg %>% left_join(maj_annot, by = "gut_cell_type")
saveRDS(deg, paste0(crohns_dir, "/deg/crohns_deg_pre-processed_revision.RDS"))

# Change the major cell type annotation for ILCs to NK instead of Unconventional T and for Plasmablast to Plasma B instead of B
# NOTE: metadata/cell.tsv is the authoritative cell map and already carries this revision in its
# revision_major_cell_type column (ILC -> NK, Plasmablast -> Plasma B), so it is read directly below
# and no revised cell map is written.
######## 

# Note the shared cell types 
tenk_maj_annot <- read.delim("metadata/cell.tsv")
common_gut_cell_types <- intersect(tenk_maj_annot$revision_major_cell_type, maj_annot$major_cell_type)
write.table(common_gut_cell_types, paste0(crohns_dir, "/deg/common_gut_cell_types_with_tenk_revision.txt"), sep = "\t", col.names = F, quote = F, row.names = F)
