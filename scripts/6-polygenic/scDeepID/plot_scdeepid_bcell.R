# Extended Data Fig 7: scDeepID cell function scores in SLE-enriched B cells.
# a) UMAP of the scDeepID embedding, coloured by five cell functions and by cell type;
# b) the scores in SLE-enriched (scDRS Monte Carlo P < 0.05) and other cells, with
#    the one-sided t-test P; c) mean score by decile of the scDRS SLE score.
# Ported from scripts/4-polygenic/scDeepID/scDRS_scDeepID_BCELL.py (main branch).
# The per-cell scores come from rule scdeepid_bcell_extract (workflow/rules/scdrs.smk);
# the P values from the pathway test behind Supp Table 10.

library(tidyverse)
library(arrow)
library(patchwork)
library(ggrastr)
library(ragg)

group_colours <- c("Non-SLE" = "lightblue", "SLE" = "lightcoral")
line_colours <- c("lightcoral", "lightskyblue", "lightgreen", "peachpuff", "plum")
pathway_label <- function(x) {
  str_remove(x, "^GOBP_") |>
    str_replace_all("_", " ") |>
    str_to_sentence() |>
    str_replace_all(c("Nf kappab" = "NF-κB", "nf kappab" = "NF-κB", "erk1 and erk2" = "ERK1 and ERK2"))
}

df_cells <- read_parquet("results/scDeepID/tenk10k_bcell.sle.ed_fig7_pathways.anon.parquet")
pathways <- setdiff(names(df_cells), c("index", "umap_1", "umap_2", "scdeepid_cell_type"))
labels <- set_names(str_wrap(pathway_label(pathways), 32), pathways)

df_cells <- df_cells |>
  inner_join(read_parquet("results/aggregate/tenk10k_phase1.scdrs.cell_score.anon.parquet.gz",
                          col_select = c("index", "sle")) |> rename(sle_score = sle), by = "index") |>
  inner_join(read_parquet("results/aggregate/tenk10k_phase1.scdrs.cell_mcp.anon.parquet.gz",
                          col_select = c("index", "sle")) |> rename(sle_mcp = sle), by = "index") |>
  mutate(sle_status = factor(if_else(sle_mcp < 0.05, "SLE", "Non-SLE"), names(group_colours)))
stopifnot(nrow(df_cells) == 30000)

df_long <- df_cells |>
  pivot_longer(all_of(pathways), names_to = "pathway", values_to = "score") |>
  mutate(pathway = factor(pathway, pathways, labels))

df_p <- read_csv("results/scDeepID/pathway_test_full.csv", show_col_types = FALSE) |>
  filter(pathway %in% pathways) |>
  mutate(pathway = factor(pathway, pathways, labels), label = sprintf("p = %.3e", p_value))

theme_panel <- theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 7))

# a) UMAP -----------------------------------------------------------------------
umap_axes <- theme(axis.text = element_blank(), axis.ticks = element_blank())

p_umap <- map(pathways, \(pw)
  df_cells |>
    arrange(.data[[pw]]) |>
    ggplot(aes(umap_1, umap_2, colour = .data[[pw]])) +
    rasterise(geom_point(size = 0.05), dpi = 300) +
    scale_colour_viridis_c(option = "magma", direction = -1, name = NULL) +
    labs(title = labels[[pw]], x = NULL, y = NULL) +
    theme_panel + umap_axes)

centres <- df_cells |>
  group_by(scdeepid_cell_type) |>
  summarise(across(c(umap_1, umap_2), median))
p_celltype <- ggplot(df_cells, aes(umap_1, umap_2, colour = scdeepid_cell_type)) +
  rasterise(geom_point(size = 0.05), dpi = 300) +
  geom_text(data = centres, aes(label = scdeepid_cell_type), colour = "black", size = 2.2) +
  scale_colour_manual(values = c(B_intermediate = "#1f77b4", B_memory = "#ff7f0e", B_naive = "#2ca02c"),
                      guide = "none") +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_panel + umap_axes

panel_a <- wrap_plots(c(p_umap, list(p_celltype)), ncol = 3)

# b) SLE-enriched and other cells ------------------------------------------------
panel_b <- ggplot(df_long, aes(sle_status, score)) +
  geom_violin(aes(fill = sle_status), colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_text(data = df_p, aes(x = 1.5, y = Inf, label = label), vjust = 1.5, size = 2.2) +
  facet_wrap(~ pathway, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = group_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Cell function score") +
  theme_panel +
  theme(strip.background = element_blank(), strip.text = element_text(size = 7, hjust = 0))

# c) by scDRS score decile --------------------------------------------------------
panel_c <- df_long |>
  group_by(pathway) |>
  mutate(decile = ntile(sle_score, 10)) |>
  group_by(pathway, decile) |>
  summarise(score = mean(score), .groups = "drop") |>
  ggplot(aes(decile, score, colour = pathway)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_x_continuous(breaks = 1:10, labels = paste0("Q", 1:10)) +
  scale_colour_manual(values = line_colours, name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(x = "scDRS Score Quantiles", y = "Mean Pathway Score") +
  theme_bw(base_size = 7) +
  theme(legend.position = "top", legend.direction = "vertical", legend.justification = "left",
        legend.text = element_text(size = 5), legend.key.size = unit(2.5, "mm"),
        legend.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank())

tag <- theme(plot.tag = element_text(face = "bold", size = 8))
p <- (wrap_elements(panel_a) + labs(tag = "a") + tag) /
  (panel_b + labs(tag = "b") + tag + inset_element(panel_c + labs(tag = "c") + tag, 0.69, 0, 1, 0.5))

ggsave("figures/fig-e7.png", p, width = 7.1, height = 9.4, device = agg_png, dpi = 300)
ggsave("figures/fig-e7.pdf", p, width = 7.1, height = 9.4, device = cairo_pdf)
