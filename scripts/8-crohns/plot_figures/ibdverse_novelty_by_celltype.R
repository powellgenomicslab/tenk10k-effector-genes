# Supplementary Figure 25: enrichment of IBDverse-nominated effector genes already found in
# TenK10K, by IBDverse cell type, against the fraction of blood-derived cells.
# Reproduces the cell-type figure from
# TenK10K_causal_inference_vs_IBDverse/scripts/enrichment.r (submodule), with
# larger text and points for publication.

library(paletteer)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
library(scales)

input.dir <- "results/pi1_mr_novelty_results"
out.dir <- "figures"
out.name <- "fig-s25"

# text/point sizes (original figure: label = 6, point = 4, strip = 7, base = 11, 9 x 6 in)
sz <- list(label = 9, point = 5, strip = 10, base = 13, width = 12, height = 9)

# min_x = pseudocount for cell type with prop_blood = 0
min_x <- -5

tenk_not_novel_rel_res_celltype <- fread(paste0(input.dir, "/tenk_not_novel_res_relative_celltype.tsv")) |> 
    filter(is.finite(OR), is.finite(prop_blood)) |> 
    mutate(prop_blood_mod = ifelse(prop_blood == 0, 10^min_x, prop_blood))
tenk_not_novel_rel_res_celltype[order(category_new, label_machine),  cell_num := seq_along(label_machine)]

lm_fit <- lm(log(OR) ~ log10(prop_blood_mod), data = tenk_not_novel_rel_res_celltype)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood_mod)", "Pr(>|t|)"], 3)
print(paste0("Slope: ", slope, ". p=", pval))

p_main <- tenk_not_novel_rel_res_celltype |> 
  ggplot(aes(x = log10(prop_blood_mod), y = log(OR))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = cell_num, color = category_new), size = sz$label/.pt, seed = 42) +
  geom_point(aes(fill = category_new), size = sz$point, shape = 21, stroke = 0.7) +
  labs(
    x = "Fraction of cells derived from blood samples in IBDverse",
    y = "Odds ratio for IBDverse-nominated \neffector gene already found in TenK10K"
  ) +
  theme_classic(base_size = sz$base) +
  scale_fill_paletteer_d("ggsci::category10_d3", name = NULL) +
  scale_color_paletteer_d("ggsci::category10_d3", name = NULL) +
  annotate("text", x = -Inf, y = Inf, label = paste0("β = ", slope, ", p = ", pval),
           hjust = -0.2, vjust = 1.2, size = sz$base / .pt) +
  scale_y_continuous(
    breaks = log(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)) +
  scale_x_continuous(
    breaks = c(min_x, -4, -3, -2, -1, log10(0.5)), 
    labels = c(0, "0.01%", "0.1%", "1%", "10%", "50%")) +
  theme(legend.position = "none")

p_index <- tenk_not_novel_rel_res_celltype |>
  arrange(category_new, cell_num) |>
  mutate(row = row_number()) |>
  ggplot(aes(x = 0, y = -row * 1.5)) +
  geom_text(aes(label = paste0(cell_num, ": ", label_new), color = category_new),
            hjust = 0, size = sz$label/.pt, lineheight = 1) +
  scale_color_paletteer_d("ggsci::category10_d3") +
  facet_wrap(~ category_new, ncol = 1, scales = "free_y", space = "free_y", strip.position = "top") +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = sz$strip, face = "bold", angle = 0, hjust = 0, margin = margin(b = 0)),
    panel.spacing.y = unit(0.25, "lines"),
    strip.clip = "off"
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(expand = expansion(add = c(1)))

plots_or_celltype <- (p_main | p_index) + plot_layout(widths = c(2.4, 1))

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out.dir, paste0(out.name, ".pdf")), plots_or_celltype,
       width = sz$width, height = sz$height, device = cairo_pdf, bg = "white")
ggsave(file.path(out.dir, paste0(out.name, ".png")), plots_or_celltype,
       width = sz$width, height = sz$height, device = ragg::agg_png, dpi = 300, bg = "white")
