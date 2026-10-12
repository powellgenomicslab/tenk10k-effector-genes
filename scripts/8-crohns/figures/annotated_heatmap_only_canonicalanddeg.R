library(tidyverse)
library(patchwork)

res <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_all.RDS") %>% 
  ungroup() %>% 
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR)) 


canon_genes <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/sig_canon_genes.RDS") 
eqtl_gen_genes <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/eqtlgen_crohns_sig_genes.RDS") 
drug_genes <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/sig_drug_pathway_or_target_genes.RDS")
ct_spec_genes <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/tenk_crohns_sig_cell_type_specific_genes.RDS")

# load deg genes lists
load("resources/crohns_case_study/crohns_relevant_gene_lists/mr_deg_gene_lists.RData") # contains the full, the concordant only, discordant only gene lists
deg_genes <- concordant_genes # here we use only concordant

opp_dir_genes <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/sig_opp_dir_genes.RDS")

#######################################################################################################################################################################################################

# prepare a df for each situation 
# df1 = canon or drug

test <- res %>% filter(probeID %in% drug_genes | probeID %in% canon_genes)
test <- test %>% complete(Gene, cell_type)
df1 <- test

df1$annot = "Canonical or Drug Pathway"

##############################################################################

test <- res %>% filter(probeID %in% deg_genes) 
# & !(probeID %in% df1$probeID)) # drop genes from previous categories

test <- test %>% complete(Gene, cell_type)

df2 = test

df2$annot = "Differentially Expressed in Crohn's Disease"

# ##############################################################################
# # df3 = opp direction
# test <- res %>% filter(Gene %in% opp_dir_genes) %>% filter(!(probeID %in% df1$probeID)& !(probeID %in% df2$probeID)) # we don't want genes from any previous categories
# 
# test <- test %>% complete(Gene, cell_type)
# df3 <- test
# 
# df3$annot = "Discordant bMR"
# # If we are not including this:
# df3 = NULL
# 
# ##############################################################################
# 
# # df4 = cell type specific 
# test <- res %>% filter(probeID %in% ct_spec_genes & !(probeID %in% eqtl_gen_genes$probeID) & !(probeID %in% df1$probeID) & !(probeID %in% df2$probeID) & !(probeID %in% df3$probeID))
# 
# # get top res from ct spec
# top_genes <- test %>% filter(sig == T) %>% slice_min(p_SMR_multi, n = 20) %>% pull(probeID) %>% unique()
# 
# test <- test %>% filter(probeID %in% top_genes)
# #test <- test %>% filter(sig == TRUE) %>% sample_n(30)
# 
# test <- test %>% complete(Gene, cell_type)
# 
# df4 = test
# df4$annot = "Top Cell Type Specific & Not in Bulk"
# 
# # If we are not including this:
# df4 = NULL
##############################################################################

#plot_data = rbind(df1, df2, df3, df4)
plot_data = rbind(df1, df2)
plot_data = plot_data %>% filter(gene_type == "protein_coding")

unique(plot_data$Gene) # check how many genes you will plot 

max_ptransform <- plot_data %>% pull(p_transform) %>% abs() %>% max(na.rm = TRUE) # for a universal one! 

(FIG1 <- ggplot(plot_data, aes(y = fct_reorder(Gene, sig, ~sum(.x, na.rm = T)), x = fct_reorder(cell_type, as.integer(major_cell_type)), fill = p_transform)) +
    facet_grid(rows = vars(annot), scales = "free", space = "free") +
    # geom_tile(data = ~filter(.x, is.na(p_transform))) +
    geom_tile(fill = "grey90", colour = "white", linewidth = 0.5) +
    geom_tile(data = ~filter(.x, !is.na(p_transform)), color = "black", linewidth = 1) +
    # geom_tile(data = ~filter(.x, !is.na(p_transform)), color = "black", linewidth = 0.5) +
    geom_point(aes(shape = "MR Gene"), size = 1, 
               data = ~filter(.x, sig))+
    paletteer::scale_fill_paletteer_c("ggthemes::Red-Blue-White Diverging",
                                      na.value = "grey90", limits = c(-max_ptransform, max_ptransform), direction = -1) +
    scale_x_discrete() +
    labs(colour = NULL, x = NULL, fill = bquote(-log[10]~italic(P) %*% "direction of effect"), size = 5, shape = NULL) +
    # theme_minimal() +
    theme(text = element_text(family = "Roboto"),
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

# ggsave("causal_inference_manuscript/figs/main_corrected/heat_facet_annot_2_categories_Roboto.png", FIG1, device = ragg::agg_png(),
#        width = 5.0, height = 16, bg = "white", scaling = 1.2, dpi = 300)


save.image("resources/crohns_case_study/figures/heatmap_objects.RData")
