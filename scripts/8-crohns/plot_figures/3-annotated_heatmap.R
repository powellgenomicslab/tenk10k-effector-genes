library(tidyverse)
library(patchwork)

tenk_raw <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_raw.RDS") %>% 
  ungroup() %>% 
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR)) 

mr_deg <- readRDS("resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS")
deg_genes <- unique(mr_deg$Gene)
#######################################################################################################################################################################################################

# Create annot col for faceting
preplot_data <- tenk_raw %>% 
  mutate(
    across(where(is.logical), function(x) replace_na(x, FALSE)),
    canon = (liuetal | sys_review | drug_target | drug_pathway_tnf | drug_pathway_il23 | drug_pathway_integrin),
    deg_matched_ct = Gene %in% deg_genes,
    annot = case_when(
      deg_matched_ct & canon ~ "DEG + Literature-Annotated",
      !deg_matched_ct & canon  ~ "Literature-Annotated",
      deg_matched_ct & !canon ~ "DEG",
      TRUE ~ "Neither")) %>% 
  filter(gene_type == "protein_coding") %>% 
  mutate(annot = factor(annot, levels = c("Literature-Annotated", "DEG", "DEG + Literature-Annotated"))) %>% 
  arrange(annot)

# pull genes 
genes_to_plot <- preplot_data %>% filter(annot != "Neither" & mr_sens_coloc == TRUE) %>% pull(Gene)

# Complete so na values for p_transform indicate missing values 
plot_data <- preplot_data %>% complete(Gene, cell_type)
 # Get rid of NA annotations 
plot_data <- plot_data %>% filter(!is.na(annot))
#plot_data <- plot_data %>% filter(!is.na(annot)) %>% filter(mr_sens == TRUE & coloc == TRUE)

max_ptransform <- plot_data %>% pull(p_transform) %>% abs() %>% max(na.rm = TRUE) # for a universal one! 

# (FIG1 <- ggplot(plot_data %>% filter(Gene %in% genes_to_plot, 
#                                      aes(y = fct_reorder(Gene, mr_sens, ~sum(.x, na.rm = T)), x = fct_reorder(cell_type, as.integer(major_cell_type)), fill = p_transform))) +
(FIG1 <- ggplot(plot_data %>% filter(Gene %in% genes_to_plot), aes(y = Gene, x = fct_reorder(cell_type, as.integer(major_cell_type)), fill = p_transform)) +
    coord_fixed() +
    facet_grid(rows = vars(annot), scales = "free", space = "free") +
    # geom_tile(data = ~filter(.x, is.na(p_transform))) +
    geom_tile(fill = "grey", colour = "white", linewidth = 0.5) +
    geom_tile(data = ~filter(.x, !is.na(p_transform)), color = "black", linewidth = 0.5) +
    # geom_tile(data = ~filter(.x, !is.na(p_transform)), color = "black", linewidth = 0.5) +
    geom_point(aes(shape = "MR + Sensitivity + Coloc"), size = 1,
               data = ~filter(.x, mr_sens_coloc)) +
    paletteer::scale_fill_paletteer_c("ggthemes::Red-Blue-White Diverging",
                                      na.value = "grey", limits = c(-max_ptransform, max_ptransform), direction = -1) +
    scale_x_discrete() +
    labs(colour = NULL, x = NULL, fill = bquote(-log[10]~italic(P) %*% "direction of effect"), size = 5, shape = NULL) +
    #theme_minimal() +
    scale_y_discrete(limits = rev) +
    # fixed legend order (otherwise it can change between runs)
    guides(fill = guide_colourbar(order = 1), shape = guide_legend(order = 2)) +
    theme(text = element_text(family = "Helvetica"),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 9, vjust = 0.5, hjust = 1, face = "italic"),
          plot.title = element_text(size = 50, face = "bold"),
          legend.position = "bottom",
          panel.grid = element_blank(),
          # legend.justification = c(0, 1),
          legend.title = element_text(size = 7)) 
  # axis.text.y = element_text(face = "bold")) +
  # labs(x = NULL, shape = NULL, fill = bquote(-log[10]~italic(P) %*% "dir"))
  # facet_grid(annot) +
  # theme(strip.placement = "inside") +
  # ggtitle(label = "Annotated Causal Genes for Crohn's Disease")# Add significance layer
)

# Figure 6b panel
dir.create("figures/panels/8-crohns", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/panels/8-crohns/fig6b_heatmap_canonical_deg_mrsenscoloc.png", FIG1, device = ragg::agg_png(),
        width = 5.0, height = 16, bg = "white", scaling = 1.0, dpi = 500)
ggsave("figures/panels/8-crohns/fig6b_heatmap_canonical_deg_mrsenscoloc.pdf", FIG1, device = cairo_pdf,
       width = 5.0, height = 16, bg = "white")

save.image("resources/crohns_case_study/figures/heatmap_objects.RData")
