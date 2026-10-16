# Figure 6d: enrichment of IBDverse-nominated effector genes already found in TenK10K,
# by IBDverse major cell category, against the fraction of cells derived from blood.
# Reproduces the "relative" category figure of
# TenK10K_causal_inference_vs_IBDverse/scripts/enrichment.r (branch new_annot), whose
# inputs are written by that repository's scripts/summarise.r. The labels are placed
# with a fixed seed and no time limit, so the panel is the same on every run.

library(tidyverse)
library(data.table)
library(ggrepel)

input.dir <- "results/pi1_mr_novelty_results"
df_prop_blood <- fread(paste0(input.dir, "/prop_blood_ibdverse.csv"))
tenk_not_novel_rel_res <- fread(paste0(input.dir, "/tenk_not_novel_res_relative.tsv")) |>
    filter(is.finite(OR)) |>
    left_join(df_prop_blood, by = c("category" = "predicted_category"))

# Are IBDverse cat genes RELATIVELY enriched for non-novelty vs other IBDverse categories?
lm_fit <- lm(log(OR) ~ log10(prop_blood), data = tenk_not_novel_rel_res)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood)", "Pr(>|t|)"], 3)
print(paste0("Slope: ", slope, ". p=", pval))
min_x <- -5

tenk_not_novel_rel_res <- tenk_not_novel_rel_res %>% mutate(blood_derived = category %in% c("B", "Plasma", "Myeloid", "T"), TRUE)

p <- ggplot(tenk_not_novel_rel_res, aes(x = log10(prop_blood), y = log(OR))) +
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf,
           fill = "#4DAF4A", alpha = 0.2) +
  annotate("rect", xmin = -2, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#377EB8", alpha = 0.2) +
  annotate("text", x = -Inf, y = -Inf, label = "Gut-specific (<1% blood-derived cells)",
           hjust = -0.1, vjust = -0.5, color = "#4DAF4A", size = 3.5, fontface = "plain") +
  annotate("text", x = Inf, y = -Inf, label = "Non-gut-specific",
           hjust = 1.1, vjust = -0.5, color = "#377EB8", size = 3.5, fontface = "plain") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(fill = blood_derived), size = 4, shape = 21, stroke = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#B03060", "FALSE" = "darkblue")) +
  geom_label_repel(aes(label = category), size = 9/.pt, force = 5, min.segment.length = 0, border.colour = NA,
    label.size = 0, max.overlaps = Inf, label.padding = 0.1, point.padding = unit(5, "lines"), fill = alpha("white", .8),
    seed = 1, max.time = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  labs(
    x = "Fraction of cells derived from blood samples in IBDverse",
    y = "Odds ratio for IBDverse-nominated \neffector gene already found in TenK10K"
  ) +
  theme_classic() +
  scale_y_continuous(
    breaks = log(c(1, 1.5, 2, 3, 4, 5)),
    labels = exp) +
  scale_x_continuous(
    breaks = c(min_x, -4, -3, -2, -1, log10(0.5)),
    labels = c(0, "0.01%", "0.1%", "1%", "10%", "50%")) +
  annotate("label", x = -Inf, y = Inf, label = paste0("β = ", slope, ", p = ", pval), hjust = -0.2, vjust = 1.2) +
  theme(legend.position = "none")

out.dir <- "figures/panels/8-crohns"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out.dir, "fig6d_ibdverse_blood_enrichment.png"), p,
       device = ragg::agg_png, scaling = 1, width = 5.5, height = 5, dpi = 500)
ggsave(file.path(out.dir, "fig6d_ibdverse_blood_enrichment.pdf"), p,
       device = cairo_pdf, width = 5.5, height = 5)
