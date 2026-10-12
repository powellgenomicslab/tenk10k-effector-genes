library(tidyverse)
library(data.table)
library(patchwork)
library(ragg)
library(scales)
library(geomtextpath)

source("scripts/0-preprocess/preprocess_results.R")

# Extended Data figure: overlap of TenK10K single-cell MR with MAGMA (a) and
# eQTLgen bulk whole blood MR (b), by number of cell types with MR associations.
# Previously panels a and d of Fig. 3 (4_pub_combined.R).

mk_plot_tally <- function(group_col, label = list(mr_only = "MR only", mr_other = "MR & MAGMA"),
                          p_title = NULL,
                          col_scheme = list(mr_only = "#A4ABB0FF", mr_other = "#4C6C94FF")){
  df_tally <- df_msmr %>%
    filter(!is.na({{group_col}})) %>%
    group_by(probeID, phenotype, {{group_col}}) %>%
    tally(name = "n_celltypes") %>%
    group_by(n_celltypes, {{group_col}}) %>%
    tally(name = "n") %>%
    ungroup() %>%
    complete(n_celltypes, {{group_col}}, fill = list(n = 0)) %>%
    mutate(sum_n = sum(n), prop = n / sum_n, .by = n_celltypes) %>%
    mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))

  ggplot(df_tally, aes(x = n_celltypes, y = n)) +
    geom_col(aes(fill = {{group_col}}), width = 1, linewidth = 0.5,
             position = "stack", color = "black") +
    scale_fill_manual(values = c(col_scheme$mr_only, col_scheme$mr_other),
                      name = NULL,
                      label = c("FALSE" = label$mr_only,
                                "TRUE" = label$mr_other),
                      guide = guide_legend(ncol = 1)) +
    theme_classic() +
    labs(y = "N gene-trait MR associations",
         x = "N cell types with MR associations",
         title = p_title
         ) +
    geom_line(aes(group = {{group_col}}, color = {{group_col}}, y = prop_scaled),
              data = ~filter(.x, !{{group_col}})) +
    geom_textline(aes(group = {{group_col}},
                      label = paste("%", ifelse({{group_col}}, label$mr_other, label$mr_only)),
                      color = {{group_col}}, y = prop_scaled),
                  data = ~filter(.x, !{{group_col}}),
                  size = 9/.pt,
                  vjust = -1, hjust = 0.4) +
    geom_point(aes(color = {{group_col}}, y = prop_scaled), data = ~filter(.x, !{{group_col}})) +
    scale_color_manual(values = c("red3", "red3"),
                       guide = "none") +
    scale_x_continuous(breaks = seq(1,28,2),
                       expand = expansion(add = 0.25)) +
    scale_y_continuous(labels = ~ifelse(.x <1e3, .x, number(.x / 1e3, 1, suffix = "k")),
                       n.breaks = 10,
                       expand = expansion(mult = c(0, 0.05)),
                       sec.axis = sec_axis(
                         ~rescale(., to = c(0,1)),
                         name = NULL,
                         labels = label_percent()),
                       guide = guide_axis(cap = "none")) +
    coord_cartesian(clip = "off") +
    theme(legend.position.inside = c(0.95,1),
          legend.position = "inside",
          plot.title = element_text(size = 10, hjust = 0, margin = margin(l = 12, b = 5.5)),
          plot.title.position = "plot",
          axis.title = element_text(size = 9),
          legend.text = element_text(margin = margin(l = 2)),
          legend.justification = c(1, 1),
          legend.key.size = unit(0.75, "lines"),
          legend.key.spacing.y = unit(0.5, "lines"),
          axis.text = element_text(size = 8),
          axis.line.y.right = element_line(color = "red3"),
          axis.text.y.right = element_text(color = "red3"),
          axis.ticks.y.right = element_line(color = "red3"),
          axis.title.y.right = element_text(color = "red3"),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
    )
}

df_msmr <- df_msmr_tenk10k |> filter(mr == TRUE)

(p_mr_gwas <- mk_plot_tally(magma_gene, p_title = "TenK10K single-cell MR vs. MAGMA"))
(p_mr_eqtlgen <- mk_plot_tally(eqtlgen_mr,
                              p_title = "TenK10K single-cell MR vs. eQTLgen bulk whole blood MR",
                        label = list(mr_only = "TenK10K only", mr_other = "TenK10K & eQTLgen"))
)

plots_overlap <- wrap_plots(
  list(wrap_elements(full = p_mr_gwas),
       wrap_elements(full = p_mr_eqtlgen)),
  ncol = 2
) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

# Extended Data Fig. 4a-b
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig-e4.pdf", plots_overlap,
       width = 8, height = 3, device = cairo_pdf)
ggsave("figures/fig-e4.png", plots_overlap,
       width = 8, height = 3, device = ragg::agg_png, dpi = 300)
