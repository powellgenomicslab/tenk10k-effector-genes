# ==========================================================================================================================================
# Fig 2: SMR overview
# ==========================================================================================================================================

suppressPackageStartupMessages(suppressWarnings({
  library(patchwork)
  library(ragg)
  library(ComplexHeatmap)
  library(tidyverse)
}))

# Panels are read from RDS files written by (rerun these first when they change):
#   a  heatmap            scripts/2-mr/trait-level-analysis/2-fig_2_heatmap.R
#   b  pleiotropy boxplots (split by 1 / 2+ cell types)  scripts/2-mr/trait-level-analysis/3-gene_features_pleiotropy.R
#   c  evidence tally     scripts/2-mr/trait-level-analysis/4-count_evidence.R
#   Extended Data Fig 1  scripts/2-mr/trait-level-analysis/3-gene_features_pleiotropy.R

# ------------------------------------------------------------------------------------------------------------------------------------------
# Panel A: heatmap spearman correllation
# ------------------------------------------------------------------------------------------------------------------------------------------

dir.create("figures/publication_pdf", recursive = TRUE, showWarnings = FALSE)

disease_heatmap <- readRDS("figures/mr_overview/figure_2_heatmap_disease_only.rds")
legend_disease <- readRDS("figures/mr_overview/figure_2_heatmap_legend_disease.rds")

# grab viewport = panel a slot (5/8 of the 6.27 in width, 4/7 of the 6.8 in height); a viewport larger than the
# slot is centred on it, which pushes the heatmap's top edge out of the slot and clips the first rows
ht_grab_width <- 6.27 * 5 / 8
ht_grab_height <- 6.8 * 4 / 7

p_heatmap_disease <- grid::grid.grabExpr(
  {
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
p_n_phenotypes_by_n_celltypes <- readRDS("figures/mr_overview/n_phenotypes_by_n_celltypes_split.rds")

#-------------------------------------------------------------------------------------------------------------------------------------------
# Panel C: Evidence tally
#-------------------------------------------------------------------------------------------------------------------------------------------

p_upset <- readRDS("figures/strict/evidence_tally_plot.rds")

# -------------------------------------------------------------------------------------------------------------------------------------------
# Combine into multipanel figure
# -------------------------------------------------------------------------------------------------------------------------------------------

design <- "
AAAAABBB
AAAAABBB
AAAAABBB
AAAAABBB
CCCCCCCC
CCCCCCCC
CCCCCCCC
"

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

ggsave("figures/publication_pdf/Fig2 - MR Summary.pdf",
       dpi = 320,
       plots, width = 6.27, height = 6.8, device = cairo_pdf)

ggsave("figures/publication_pdf/Fig2 - MR Summary.png",
       dpi = 400,
       plots, width = 6.27, height = 6.8, device = ragg::agg_png)

# -------------------------------------------------------------------------------------------------------------------------------------------
# Extended Data Fig 1: scatter plot of pleiotropic, cell type specific genes
# -------------------------------------------------------------------------------------------------------------------------------------------

p_scatter_pleiotropic <- readRDS("figures/mr_overview/fig_1_pleiotropic_genes_scatter.rds")

ggsave("figures/publication_pdf/Extended Data Fig1 - Pleiotropic genes.pdf",
       dpi = 320,
       p_scatter_pleiotropic, width = 9, height = 3.5, device = cairo_pdf)

ggsave("figures/publication_pdf/Extended Data Fig1 - Pleiotropic genes.png",
       dpi = 400,
       p_scatter_pleiotropic, width = 9, height = 3.5, device = ragg::agg_png)
