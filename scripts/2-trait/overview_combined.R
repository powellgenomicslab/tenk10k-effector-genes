# ==========================================================================================================================================
# Fig 2: SMR overview
# ==========================================================================================================================================

suppressPackageStartupMessages(suppressWarnings({
  library(patchwork)
  library(ragg)
  library(scales)
  library(paletteer)
  library(ggnewscale)
  library(geomtextpath)
  library(ComplexHeatmap)
  library(data.table)
  library(tidyverse)
}))

# Panels are read from RDS files written by (rerun these first when they change):
#   a  heatmap            scripts/2-trait/fig_2_heatmap.R
#   b  pleiotropy boxplots (split by 1 / 2+ cell types)  scripts/2-trait/gene_features_pleiotropy.R
#   c  evidence tally     scripts/3-pleiotropy/count_evidence.R
#   Extended Data Fig 1  scripts/2-trait/gene_features_pleiotropy.R

# ------------------------------------------------------------------------------------------------------------------------------------------
# Panel A: heatmap spearman correllation
# ------------------------------------------------------------------------------------------------------------------------------------------

disease_heatmap <- readRDS("figures/panels/2-trait/figure_2_heatmap_disease_only.rds")
legend_disease <- readRDS("figures/panels/2-trait/figure_2_heatmap_legend_disease.rds")

PLOT_SIZE <- 0.7
cell_size_cm <- unit(0.2 * PLOT_SIZE, "cm")
# grab viewport = panel a slot (5/8 of the 6.27 in width, 4/7 of the 6.8 in height); a viewport larger than the
# slot is centred on it, which pushes the heatmap's top edge out of the slot and clips the first rows
ht_grab_width <- 6.27 * 5 / 8
ht_grab_height <- 6.8 * 4 / 7

p_heatmap_disease <- grid::grid.grabExpr(
  {
    ht_opt$DENDROGRAM_PADDING = unit(0.7 * PLOT_SIZE, "cm") # increase padding to fit in long labels
    # no vertical padding (bottom, left, top, right): the layout must fit the slot or ComplexHeatmap centres it and clips the top rows
    draw(disease_heatmap, show_heatmap_legend = FALSE, padding = unit(c(0, 2, 0, 2), "mm"))
    # draw legend in the top right corner, adjust the x and y values to position it correctly
    decorate_heatmap_body("Spearman", {
    draw(legend_disease, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
    })
  },
  width = ht_grab_width,
  height = ht_grab_height
)

#-------------------------------------------------------------------------------------------------------------------------------------------
# Panel B: Pleitropic genes box plot
#-------------------------------------------------------------------------------------------------------------------------------------------

# split version (genes in 1 cell type by cell type, over genes in 2+ cell types); text sizes are set in the panel script
p_n_phenotypes_by_n_celltypes <- readRDS("figures/panels/2-trait/n_phenotypes_by_n_celltypes_split.rds")


#-------------------------------------------------------------------------------------------------------------------------------------------
# Panel C: Evidence tally
#-------------------------------------------------------------------------------------------------------------------------------------------

p_upset <- readRDS("figures/panels/3-pleiotropy/evidence_tally_plot.rds")

# -------------------------------------------------------------------------------------------------------------------------------------------
# Combine into multipanel figure
# -------------------------------------------------------------------------------------------------------------------------------------------

# layout <- c(
#   area(t = 1, l = 1, b = 4, r = 4),
#   area(t = 1, l = 4, b = 4, r = 6),
#   area(t = 4, l = 1, b = 6, r = 3),
#   area(t = 4, l = 3, b = 6, r = 6)
# )

design <- "
AAAAABBB
AAAAABBB
AAAAABBB
AAAAABBB
CCCCCCCC
CCCCCCCC
CCCCCCCC
"


# design <- "
# AABBD
# AABBD
# CCC#D
# "

plots <- list(
  A = wrap_elements(full = p_heatmap_disease),
  B = wrap_elements(full = p_n_phenotypes_by_n_celltypes),
  C = wrap_elements(full = p_upset)
) %>%
  wrap_plots(design = design, guides = "keep") +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin())) &
  theme(
    plot.tag = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    plot.tag.position = c(0.01, 0.99)
  )

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig-m2.pdf",
       dpi = 320,
       plots, width = 6.27, height = 6.8, device = cairo_pdf)

ggsave("figures/fig-m2.png",
       dpi = 400,
       plots, width = 6.27, height = 6.8, device = ragg::agg_png)

# -------------------------------------------------------------------------------------------------------------------------------------------
# Extended Data Fig 1: scatter plot of pleiotropic, cell type specific genes
# -------------------------------------------------------------------------------------------------------------------------------------------

p_scatter_pleiotropic <- readRDS("figures/panels/2-trait/fig_1_pleiotropic_genes_scatter.rds")

ggsave("figures/fig-e1.pdf",
       dpi = 320,
       p_scatter_pleiotropic, width = 9, height = 3.5, device = cairo_pdf)

ggsave("figures/fig-e1.png",
       dpi = 400,
       p_scatter_pleiotropic, width = 9, height = 3.5, device = ragg::agg_png)


# ------------------------------------------------------------------------------------------------------------------------------------------
# Reported STATS
# ------------------------------------------------------------------------------------------------------------------------------------------

# # get overall numbers for gwas / eqtlgen annotation
# (df_tally_gwas <- df_msmr %>% 
#   group_by(probeID, phenotype, magma_gene) %>%
#   tally(name = "n_celltypes") %>%
#   group_by(magma_gene) %>%
#   tally(name = "n") %>%
#   ungroup() %>% 
#   mutate(sum_n = sum(n), prop = n / sum_n)
# )

# (df_tally_eqtlgen <- df_msmr %>% 
#     group_by(probeID, phenotype, eqtlgen_mr) %>%
#     tally(name = "n_celltypes") %>%
#     group_by(eqtlgen_mr) %>%
#     tally(name = "n") %>%
#     ungroup() %>% 
#     mutate(sum_n = sum(n), prop = n / sum_n)
# )

# # bins trait based on n association
# df_trait_category <- df_smr_pheno_cell %>%
#   filter(n_mr > 1) %>%
#   summarise(n_mr_trait = mean(n_mr), .by = c(pheno_label)) %>%
#   mutate(ntile_group = ntile(n_mr_trait, 4)) %>%
#   mutate(min_n_mr_trait = min(n_mr_trait),
#          max_n_mr_trait = max(n_mr_trait),
#          .by = ntile_group)

# trait_cat_col <- paletteer_d("MetBrewer::Hiroshige", 10) %>%
#   set_names(cat_order)

# mk_plot <- function(df_plot, n_breaks = 10) {
#   df_plot <- df_plot %>% 
#     mutate(fct_pheno = fct_reorder(pheno_label, n_mr, sum))
  
#   (p_grid <- df_plot %>% 
#     ggplot(aes(y = cell_type, x = fct_pheno)) +
#     theme_minimal() +
#     # coord_fixed() +
#     facet_grid(rows = vars(major_cell_type), cols = vars(pheno_cat),
#                scales = "free", space = "free") +
#     geom_tile(color = "black", aes(fill = n_mr), linewidth = 0.5) +
#     scale_fill_viridis_b(
#       n.breaks = 10, direction = -1, option = "magma",
#       name = "Number of MR evidence"
#     ) +
#     scale_y_discrete(expand = expansion(0.0)) +
#     scale_x_discrete(expand = expansion(add = 0.5)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7, vjust = 0.5),
#           strip.clip = "off",
#           strip.text = element_blank(),
#           plot.margin = margin())
#   )
  
#   (p_count <- df_plot %>% 
#       .[,list(n_mr = sum(n_mr, na.rm = TRUE)), 
#          by = .(fct_pheno, pheno_cat, supercategory)] %>% 
#       ggplot(aes(y = n_mr, x = fct_pheno, fill = pheno_cat)) +
#       theme_void() +
#       coord_cartesian(clip = "off") +
#       facet_grid(cols = vars(pheno_cat),
#                  scales = "free", space = "free") +
#       scale_fill_manual(values = trait_cat_col, name = NULL) +
#       geom_segment(
#         aes(xend = fct_pheno, yend = 0),
#         linewidth = 0.5, color = "gray30", alpha = 1
#       ) +
#       geom_point(size = 3, shape = "circle filled") +
#       layer_scales(p_grid)$x +
#       theme(strip.clip = "off")
#     )
#   ((p_count / p_grid) &
#     # scale_x_discrete(expand = expansion(0.0)) &
#     theme(panel.spacing = unit(0.5, "lines"))) +
#     plot_layout(ncol = 1, heights = c(0.2, 1))
  
# }

# setDT(df_smr_pheno_cell)
# df_smr_pheno_cell %>% 
#   inner_join(df_trait_category) %>%
#   # mutate(n_mr = na_if(n_mr, 0)) %>% 
#   filter(ntile_group == 3) %>%
#   ggplot(aes(y = cell_type, x = pheno_label)) +
#   theme_minimal() +
#   coord_fixed() +
#   # geom_blank() +
#   facet_grid(rows = vars(major_cell_type), cols = vars(pheno_cat),
#              scales = "free", space = "free") +
#   geom_tile(color = "white", aes(fill = n_mr)) +
#   scale_fill_viridis_b(
#     # colours = paletteer_d("MoMAColors::Exter", n = 10),
#     n.breaks = 10, direction = -1, option = "magma",
#     name = "Number of MR evidence"
#   ) +
#   # scale_fill_paletteer_d(
#   #   palette = "MoMAColors::Exter",
#   #   name = "Number of MR evidence"
#   # ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7, vjust = 0.5),
#         panel.border = element_rect(color = "black", fill = NA),
#         strip.clip = "off")


# ------------------------------------------------------------------------------------------------------------------------------------------
# panel a: tally by N celltypes and eQTLgen MR
# ------------------------------------------------------------------------------------------------------------------------------------------

#   # df_tally <- df_msmr %>% 
#   #   group_by(probeID, phenotype, eqtlgen_mr) %>%
#   #   tally(name = "n_celltypes") %>%
#   #   group_by(n_celltypes, eqtlgen_mr) %>%
#   #   tally(name = "n") %>% 
#   #   ungroup() %>% 
#   #   mutate(sum_n = sum(n), prop = n / sum_n, .by = n_celltypes) %>% 
#   #   mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))

#   # df_msmr$eqtlgen_mr %>% is.na() %>% table()

# mk_plot_tally <- function(group_col, label = list(mr_only = "MR only", mr_other = "MR & GWAS"),
#                           p_title = NULL,
#                           col_scheme = list(mr_only = "#A4ABB0FF", mr_other = "#4C6C94FF")){
#   df_tally <- df_msmr %>% 
#     group_by(probeID, phenotype, {{group_col}}) %>%
#     tally(name = "n_celltypes") %>%
#     group_by(n_celltypes, {{group_col}}) %>%
#     tally(name = "n") %>% 
#     ungroup() %>% 
#     mutate(sum_n = sum(n), prop = n / sum_n, .by = n_celltypes) %>% 
#     mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))
  
#   ggplot(df_tally, aes(x = n_celltypes, y = n)) +
#     geom_col(aes(fill = {{group_col}}), width = 1, linewidth = 0.5,
#              position = "stack", color = "black") +
#     scale_fill_manual(values = c(col_scheme$mr_only, col_scheme$mr_other),
#                       name = NULL,
#                       label = c("FALSE" = label$mr_only,
#                                 "TRUE" = label$mr_other),
#                       guide = guide_legend(ncol = 1)) +
#     theme_classic() +
#     labs(y = "N gene-trait MR associations",
#          x = "N cell types with MR associations",
#          title = p_title
#          ) +
#     geom_line(aes(group = {{group_col}}, color = {{group_col}}, y = prop_scaled),
#               data = ~filter(.x, !{{group_col}})) +
#     geom_textline(aes(group = {{group_col}},
#                       label = paste("%", ifelse({{group_col}}, label$mr_other, label$mr_only)),
#                       color = {{group_col}}, y = prop_scaled),
#                   data = ~filter(.x, !{{group_col}}),
#                   size = 9/.pt,
#                   vjust = -1, hjust = 0.1) +
#     geom_point(aes(color = {{group_col}}, y = prop_scaled), data = ~filter(.x, !{{group_col}})) +
#     scale_color_manual(values = c("red3", "red3"),
#                        guide = "none") +
#     scale_x_continuous(breaks = seq(1,28,4),
#                        expand = expansion(add = 0.25)) +
#     scale_y_continuous(labels = ~ifelse(.x <1e3, .x, number(.x / 1e3, 1, suffix = "k")),
#                        n.breaks = 10,
#                        expand = expansion(mult = c(0, 0.05)),
#                        sec.axis = sec_axis(
#                          ~rescale(., to = c(0,1)),
#                          name = NULL,
#                          labels = label_percent()),
#                        guide = guide_axis(cap = "none")) +
#     coord_cartesian(clip = "off") +
#     theme(legend.position.inside = c(0.95,1),
#           legend.position = "inside",
#           plot.title = element_text(size = 9, hjust = 0),
#           axis.title = element_text(size = 9),
#           legend.text = element_text(size = 9, margin = margin(l = 2)),
#           legend.justification = c(1, 1),
#           legend.key.size = unit(0.75, "lines"),
#           legend.key.spacing.y = unit(0.5, "lines"),
#           axis.text = element_text(size = 8),
#           # plot.margin = margin(t = 2, unit = "lines"),
#           axis.line.y.right = element_line(color = "red3"),
#           axis.text.y.right = element_text(color = "red3"),
#           axis.ticks.y.right = element_line(color = "red3"),
#           axis.title.y.right = element_text(color = "red3"),
#           panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
#     )
# }

# (p_mr_gwas <- mk_plot_tally(magma_gene, p_title = "TenK10K single-cell MR vs. GWAS"))
# (p_mr_eqtlgen <- mk_plot_tally(eqtlgen_mr, 
#                               p_title = "TenK10K single-cell MR vs. eQTLgen bulk whole blood MR",
#                         label = list(mr_only = "TenK10K only", mr_other = "TenK10K & eQTLgen"))
# )

# # individual figure gwas and eqtlgen comparison
# ggsave("figures/mr_overview/mr_gwas_comparison.png",
#        dpi = 320,
#        p_mr_gwas, width = 6, height = 4, device = ragg::agg_png, scaling = 1.2)
# ggsave("figures/mr_overview/mr_eqtlgen_comparison.png",
#        dpi = 320,
#        p_mr_eqtlgen, width = 6, height = 4, device = ragg::agg_png, scaling = 1.2)


# trait_cat_col <- paletteer_d("ggthemes::Tableau_10", 10) %>%
#   set_names(cat_order)
