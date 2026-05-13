library(tidyverse)

deg <- readRDS("resources/crohns_case_study/deg/scIBD.deg_by_major_cluster.rds")

deg <- dplyr::bind_rows(deg, .id = 'major_cell_type')

# unique(scibd[scibd$major_cell_type == "Myeloid",]$cluster)
# [1] "Non-classical monocyte" "Classical monocyte"     "Inflammatory monocyte"  "APOE+ macrophage"       "LYVE1+ macrophage"      "AREG+ macrophage"      
# [7] "Cycling macrophage"     "cDC1"                   "cDC2"                   "pDC"                    "LAMP3+ DC"              "Mast"                  
# [13] "Megakaryocyte" 

deg <- deg %>% 
  mutate(major_cell_type = major_cell_type %>% 
           dplyr::replace_when(str_detect(cluster, "monocyte") ~ "Monocyte", 
                               str_detect(cluster, "macrophage") ~ "macrophage", 
                               cluster %in% c("cDC1", "cDC2", "pDC") ~ "Dendritic",
                               cluster %in% c("LAMP3+ DC", "Mast", "Megakaryocyte") ~ "Other",
                               cluster %in% c("cDC1", "cDC2", "pDC") ~ "Dendritic",
                               major_cell_type == "CD4T" ~ "CD4 T",
                               major_cell_type == "CD8T" ~ "CD8 T",
                               major_cell_type == "B_Plasma" ~ "Plasma B")) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(deg_sig = TRUE) %>% 
  rename(Gene = gene) %>% 
  mutate(new_log2FC = log2(avg_log2FC)) %>% 
  group_by(Gene, major_cell_type) %>% 
  slice_min(p_val_adj, n = 1) %>% 
  distinct() %>% 
  mutate(p_transform = -log10(p_val)*sign(new_log2FC)) 

# tenk10k res
# update MR major cell type category 
mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% 
  #filter(mr_sens_coloc) %>% 
  mutate(major_cell_type = as.character(major_cell_type) %>% dplyr::replace_when(cell_type == "ILC" ~ "NK", cell_type == "Plasmablast" ~ "Plasma B")) %>% 
  group_by(Gene, major_cell_type) %>% 
  slice_min(p_SMR_multi, n = 1) %>% 
  ungroup() %>% 
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR))

print("Number of significant DEGs: ")
length(unique(deg$Gene))

# filter deg genes
# deg_filtered <- deg %>% filter(Gene %in% mr$Gene) 
# length(unique(deg_filtered$Gene))
# 135 
print("Number of intersecting MR + sig DEGs (probeID and Gene: ")
#deg %>% filter(probeID %in% mr$probeID) %>% pull(probeID) %>% unique() %>% length() 
deg %>% filter(Gene %in% mr$Gene) %>% pull(Gene) %>% unique() %>% length() 

# Extra probe IDs, probs non-coding 
############################################################################################################################################################
library(ggrepel) 
print("Inner join MR and DEG results by gene and major cell type")

# Merged the MR and DEG dataframes on gene and major cell type
mr_deg <- inner_join(mr, deg, by = c("Gene", "major_cell_type"), suffix = c("_mr", "_deg")) %>% 
  mutate(deg_maj_ct = TRUE) %>% 
  mutate(mr_deg_sign_concord = (sign(b_SMR) == sign(new_log2FC))) %>% 
  #mutate(mr_deg_sign_concord = (sign(b_SMR) == 1 & avg_log2FC > 1) | (sign(b_SMR) == -1 & avg_log2FC < 1)) %>% 
  group_by(Gene) %>% 
  mutate(mr_deg_sign_ct_het = ifelse(n_distinct(mr_deg_sign_concord) == 1, FALSE, TRUE)) %>% 
  ungroup() %>% 
  select(Gene, major_cell_type, cell_type, deg_sig, deg_maj_ct, mr_deg_sign_concord, mr_deg_sign_ct_het, everything()) %>% 
  distinct()

print(mr_deg)

# Print statistics for writing 

print("Number of unique genes MR significant results: ")
length(unique(mr$Gene))

print("Number of unique genes deg filtered to significant results: ")
length(unique(deg$Gene))

print("Number of rows in mr_deg: ")
n_distinct(mr_deg) 

print("Number of unique genes in mr_deg: ")
length(unique(mr_deg$Gene))

print("Number of concordant gene-pbmc cell type-gut cell type combos: ")
mr_deg %>% group_by(mr_deg_sign_concord) %>% tally()

print("Among the concordant mr_deg results at the major cell type level, how many were discordant between cell types?")
mr_deg %>% filter(mr_deg_sign_concord == TRUE) %>% group_by(mr_deg_sign_ct_het) %>% tally()




darker_colors <- c("#C78E00", "#674FA3", "#1976D2", "#3E2723", "#D84315", "#1B5E20", "#FBC02D", "#4A001F") 

mr_scale <- 70
deg_scale <- 70

scatterplot_all_in_one <- ggplot(mr_deg, aes(x = p_transform_mr, y = new_log2FC, colour = major_cell_type, label = Gene)) +
  #scatterplot <- ggplot(mr_deg, aes(x = b_SMR, y = `Discrete DE coefficients`, colour = major_cell_type, label = Gene)) +
  # facet_wrap(~major_cell_type) +
  scale_colour_manual(values = darker_colors) +
  geom_point(alpha = 0.5) + 
  coord_fixed() +
  #scale_y_continuous(limits = c(-deg_scale, deg_scale), n.breaks = 10) +
  #scale_x_continuous(limits = c(-mr_scale, mr_scale), n.breaks = 10) +
  geom_text_repel(key_glyph = "point", size = 4, fontface = "italic", max.overlaps = Inf) + 
  theme_minimal() +
  geom_hline(yintercept = 0, linewidth = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.2, alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  labs(x = bquote(-log[10]~italic(p) %*% "MR direction of effect"), 
       y = "DEG log2 FC",
       colour = "Major Cell Type") + 
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 12, face = "bold")) + 
  ggtitle("Comparison of the Significance and Direction of Effect between TenK10K MR and DEG (scIBD)")

scatterplot_all_in_one

scatterplot <- ggplot(mr_deg, aes(x = p_transform_mr, y = new_log2FC, colour = major_cell_type, label = Gene)) +
  #scatterplot <- ggplot(mr_deg, aes(x = b_SMR, y = `Discrete DE coefficients`, colour = major_cell_type, label = Gene)) +
  facet_wrap(~major_cell_type, nrow = 4) +
  scale_colour_manual(values = darker_colors) +
  geom_point(alpha = 0.5) + 
  coord_fixed() +
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(n.breaks = 5) +
  geom_text_repel(size = 5, fontface = "italic", max.overlaps = Inf) + 
  theme_minimal() +
  geom_hline(yintercept = 0, linewidth = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.2, alpha = 0.5) +
  # annotate("text", x = -12, y = 12, label = "MR- DEG+", size = 3) +
  # annotate("text",x = 10, y = -12, label = "MR+ DEG-", size = 3) +
  # annotate("text", x = 10, y = 12, label = "MR+ DEG+", size = 3) +
  # annotate("text", x = -12, y = -12, label = "MR- DEG-", size = 3) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  labs(x = bquote(-log[10]~italic(p) %*% "MR direction of effect"), 
       y = bquote(-log[10]~italic(p) %*% "DEG direction of effect")) + 
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold")) +
  ggtitle("Comparison of the Significance and Direction of Effect between TenK10K MR and DEG")

scatterplot

ggsave("resources/crohns_case_study/figures/mr_deg_facet_annot_scIBD.png", scatterplot, device = ragg::agg_png(),
       width = 5.0, height = 8, bg = "white", scaling = 0.5, dpi = 300)

ggsave("resources/crohns_case_study/figures/mr_deg_scIBD.png", scatterplot_all_in_one, device = ragg::agg_png(),
       width = 5.0, height = 4.5, bg = "white", scaling = 0.7, dpi = 300)


