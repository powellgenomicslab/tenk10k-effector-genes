library(data.table)
library(tidyverse)
library(arrow)
library(ggrepel)
library(readxl)
library(RColorBrewer)

df_msmr_strict <- read_parquet("/g/data/fy54/analysis/tenk10k-causal/results/preprocessed/tenk10k_phase1.v2.parquet.gz") %>%
    setDT()

# Match trait + cell-type ordering used in other plots (see scripts/preprocess.R).
pheno_order_strict <- df_msmr_strict[, .N, by = pheno_label][order(-N), as.character(pheno_label)]
df_cell_map <- fread("resources/metadata/cell_map.tsv")
cell_type_order <- df_cell_map$cell_type

# get just genes in the MYRF / FADS locus
fads_genes <- c("FADS1", "FADS2", "FADS3", "MYRF", "TMEM258", "FEN1", "FADS1-AS1", "FADS2-AS1")
fads_results <- df_msmr_strict %>% filter(Gene %in% fads_genes)

# ------------------------------------------------------------------------------------------------
# Pleiotropic scatter plot for FADS locus genes
# (adapted from scripts/revision/gene_features_pleiotropy.R)
# ------------------------------------------------------------------------------------------------

# Keep only Gene x cell_type combinations with at least one MR-significant phenotype
sig_combos <- fads_results[mr == TRUE, .N, by = .(Gene, cell_type)][, .(Gene, cell_type)]

plot_data <- fads_results %>%
    semi_join(sig_combos, by = c("Gene", "cell_type")) %>%
    mutate(signed_neg_log10_p = -log10(p_SMR_multi) * sign(b_SMR)) %>%
    group_by(Gene, cell_type) %>%
    mutate(x = frank(-signed_neg_log10_p, ties.method = "average")) %>%
    ungroup() %>%
    mutate(
      Significant = if_else(is.na(max_evidence), "Not significant", as.character(max_evidence)),
      `Phenotype category` = pheno_cat,
    )

label_data <- plot_data %>%
  filter(mr & supercategory == "disease") %>%
  mutate(mid_x = max(x) / 2, .by = c(Gene, cell_type))

cat_order <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx",
  sheet = "trait_category_order"
) %>%
  pull(cat_order)

trait_cat_col <- brewer.pal(max(3, min(length(cat_order), 12)), "Paired")[seq_along(cat_order)] %>%
  set_names(cat_order)
disease_category_cols <- trait_cat_col[names(trait_cat_col) %in% unique(plot_data$`Phenotype category`)]

# Inside-panel strip text (top-right corner): "<cell_type>: <Gene>"
facet_label_data <- plot_data %>%
  distinct(Gene, cell_type) %>%
  mutate(label = paste0("\"", cell_type, ":\" ~ italic(\"", Gene, "\")"))

# Layout: ~4 panels per row works for ~12-30 facets
n_facets <- nrow(distinct(plot_data, Gene, cell_type))
n_col <- min(4, n_facets)
n_row <- ceiling(n_facets / n_col)

p_scatter_fads <- plot_data %>%
    ggplot(aes(shape = Significant, y = signed_neg_log10_p, x = x,
               color = `Phenotype category`, fill = `Phenotype category`)) +
    facet_wrap(vars(Gene, cell_type), scales = "free", ncol = n_col) +
    geom_point(data = plot_data %>% filter(!mr),
      color = "lightgrey", fill = "lightgrey", size = rel(2),
      show.legend = c(shape = TRUE, color = FALSE, fill = FALSE)) +
    geom_point(data = plot_data %>% filter(mr), size = rel(2)) +
    geom_text_repel(data = label_data,
      aes(label = pheno_label),
      nudge_x = label_data$mid_x - label_data$x,
      size = rel(3),
      direction = "y",
      min.segment.length = 0,
      hjust = 0.5,
      segment.color = "grey50",
      point.padding = 0.5,
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.size = 0.1,
      force = 4,
      show.legend = FALSE
    ) +
    geom_text(data = facet_label_data,
      aes(x = Inf, y = Inf, label = label),
      parse = TRUE, inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.3,
      size = 8 * 1.3 / .pt, family = "Helvetica") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    guides(
      y = guide_axis(cap = "none"),
      colour = guide_legend(ncol = 1),
      fill = guide_legend(ncol = 1),
      shape = guide_legend(ncol = 1, override.aes = list(fill = "black"))
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 8, base_family = "Helvetica") +
    labs(y = expression(sign(beta) %*% -log[10](p[SMR_multi]))) +
    scale_color_manual(values = disease_category_cols) +
    scale_fill_manual(values = disease_category_cols) +
    scale_shape_manual(
      values = c(
        "Not significant" = 3,
        "mr" = 1,
        "mr_coloc" = 2,
        "mr_sens" = 21,
        "mr_sens_coloc" = 24)) +
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "vertical",
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        strip.clip = "off")

ggsave("fig_pub/fads_locus_pleiotropic_scatter.png",
  p_scatter_fads,
  width = 5 * n_col, height = 4 * n_row,
  limitsize = FALSE)

ggsave("fig_pub/fads_locus_pleiotropic_scatter.pdf",
  p_scatter_fads,
  width = 5 * n_col, height = 4 * n_row,
  limitsize = FALSE)

# ------------------------------------------------------------------------------------------------
# Heatmap: cell_type x pheno_label, faceted by Gene, fill = sign(b_SMR) * -log10(p_SMR_multi)
# ------------------------------------------------------------------------------------------------

# Restrict to traits and cell types with at least one MR-significant FADS-gene tile
# so the panels are legible.
sig_traits     <- fads_results[mr == TRUE, unique(pheno_label)]
sig_cell_types <- fads_results[mr == TRUE, unique(cell_type)]

heatmap_data <- fads_results %>%
  filter(supercategory == "disease",
         pheno_label %in% sig_traits, cell_type %in% sig_cell_types) %>%
  mutate(
    signed_neg_log10_p = -log10(p_SMR_multi) * sign(b_SMR),
    sig                = mr & !is.na(mr),
    full_evidence      = !is.na(max_evidence) & max_evidence == "mr_sens_coloc",
    pheno_label        = factor(as.character(pheno_label), levels = pheno_order_strict),
    cell_type          = factor(as.character(cell_type),  levels = cell_type_order)
  )

max_p <- max(abs(heatmap_data$signed_neg_log10_p), na.rm = TRUE)

# Custom legend key glyphs: bordered square for MR, bordered square + dot for full evidence.
draw_key_mr <- function(data, params, size) {
  grid::rectGrob(width = unit(0.85, "npc"), height = unit(0.85, "npc"),
                 gp = grid::gpar(col = "black", fill = "grey90", lwd = 1.5))
}
draw_key_mr_sens_coloc <- function(data, params, size) {
  grid::grobTree(
    grid::rectGrob(width = unit(0.85, "npc"), height = unit(0.85, "npc"),
                   gp = grid::gpar(col = "black", fill = "grey90", lwd = 1.5)),
    grid::pointsGrob(0.5, 0.5, pch = 16,
                     gp = grid::gpar(col = "black"),
                     size = unit(0.25, "char"))
  )
}

heatmap_data_byct <- heatmap_data %>%
  mutate(Gene = factor(Gene, levels = fads_genes))

p_heatmap_fads_byct <- ggplot(heatmap_data_byct,
    aes(x = Gene, y = pheno_label, fill = signed_neg_log10_p)) +
  facet_grid(cols = vars(cell_type), scales = "free", space = "free") +
  coord_fixed(ratio = 1) +
  geom_tile(fill = "grey95", colour = "white", linewidth = 0.25) +
  geom_tile(linewidth = 0.25, colour = "grey80",
            data = ~filter(.x, !is.na(signed_neg_log10_p), !sig)) +
  geom_tile(aes(colour = "MR"), linewidth = 0.5,
            data = ~filter(.x, !is.na(signed_neg_log10_p), sig, !full_evidence),
            key_glyph = draw_key_mr) +
  geom_tile(aes(colour = "MR + SENS + COLOC"), linewidth = 0.5,
            data = ~filter(.x, full_evidence),
            key_glyph = draw_key_mr_sens_coloc) +
  geom_point(data = ~filter(.x, full_evidence),
             size = 0.6, colour = "black", show.legend = FALSE) +
  scale_colour_manual(
    values = c("MR" = "black", "MR + SENS + COLOC" = "black"),
    name   = NULL,
    guide  = guide_legend(order = 2, override.aes = list(fill = NA))) +
  scale_fill_distiller(
    palette   = "RdBu",
    direction = -1,
    na.value  = "grey90",
    limits    = c(-max_p, max_p),
    guide     = guide_colorbar(
      order          = 1,
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(7.5, "lines"),
      barheight      = unit(0.6, "lines"))) +
  labs(x = NULL, y = NULL,
       fill = expression(-log[10] ~ italic(P) %*% "direction of effect")) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(
    axis.text.x          = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                        face = "italic", size = 7),
    axis.text.y          = element_text(size = 7),
    strip.text.x         = element_text(size = 9, angle = 90,
                                        hjust = 0, vjust = 0.5),
    strip.background     = element_blank(),
    panel.background     = element_rect(fill = "grey95", colour = NA),
    panel.grid           = element_blank(),
    panel.spacing        = unit(0.4, "lines"),
    legend.position      = "bottom",
    legend.box           = "horizontal",
    legend.text          = element_text(size = 7),
    legend.title         = element_text(size = 7))

n_genes_hm_byct <- n_distinct(heatmap_data_byct$Gene)
w_hm_byct <- 13
h_hm_byct <- 5

ggsave("fig_pub/fads_locus_heatmap_by_celltype.png",
  p_heatmap_fads_byct, width = w_hm_byct, height = h_hm_byct,
  limitsize = FALSE, bg = "white")
ggsave("fig_pub/fads_locus_heatmap_by_celltype.pdf",
  p_heatmap_fads_byct, width = w_hm_byct, height = h_hm_byct, limitsize = FALSE)
