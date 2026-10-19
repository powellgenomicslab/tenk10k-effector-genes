# Supp Fig 15: mean scDeepID cell function scores in SLE-enriched (scDRS Monte Carlo
# P < 0.05) and other B cells, for the top 30 cell functions (FDR < 0.05, ranked by FDR
# of the one-sided t-test) in each of 7 GO biological-process categories.

library(tidyverse)
library(patchwork)
library(ragg)

group_colours <- c("Non-SLE" = "lightblue", "SLE" = "lightcoral")

df_test <- read_csv("results/scDeepID/pathway_test_full.csv", show_col_types = FALSE)
categories <- str_subset(names(df_test), "^cat_")

plot_category <- function(cat) {
  df <- df_test |>
    filter(.data[[cat]], fdr < 0.05) |>
    arrange(fdr, p_value) |>
    head(30) |>
    mutate(label = fct_rev(fct_inorder(str_to_sentence(str_replace_all(str_remove(pathway, "^GOBP_"), "_", " ")))))
  df |>
    pivot_longer(c(non_sle_mean, sle_mean), names_to = "group", values_to = "mean") |>
    mutate(group = factor(if_else(group == "sle_mean", "SLE", "Non-SLE"), names(group_colours))) |>
    ggplot(aes(mean, label)) +
    geom_line(aes(group = label), colour = "grey50") +
    geom_point(aes(colour = group), size = 1.5) +
    scale_colour_manual(values = group_colours, name = NULL) +
    labs(title = sprintf("%s (n=%d)", str_replace_all(str_remove(cat, "^cat_"), "_", " "), nrow(df)),
         x = NULL, y = NULL) +
    theme_bw(base_size = 6) +
    theme(plot.title = element_text(size = 6), plot.title.position = "plot", panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())
}

p <- wrap_plots(map(categories, plot_category), ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = "Mean cell function score; cell functions ranked by -log10(P) from the top") &
  theme(legend.position = "bottom", plot.caption = element_text(size = 7, hjust = 0.5))

ggsave("figures/fig-s15.png", p, width = 7.5, height = 12.5, device = agg_png, dpi = 300)
ggsave("figures/fig-s15.pdf", p, width = 7.5, height = 12.5, device = cairo_pdf)
