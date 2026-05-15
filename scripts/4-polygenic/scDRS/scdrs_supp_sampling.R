# Compare scDRS Crohn's disease results between sampling strategies

library(data.table)
library(tidyverse)
library(scales)
library(patchwork)
library(arrow)
library(ragg)
library(paletteer)
library(ggforce)
library(ggridges)
library(ggrepel)

df_cell_map <- fread("metadata/cell.tsv") %>%
    mutate(cell_type = factor(cell_type, unique(cell_type)))

df_stats <- fread("results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv")
df_top <- fread("results/aggregate/tenk10k_phase1.scdrs.cell_type_top.tsv")
df_cells <- read_parquet("results/aggregate/tenk10k_phase1/scdrs.cell_score.tsv.parquet.gz")
df_mcp <- read_parquet("results/aggregate/tenk10k_phase1/scdrs.cell_mcp.tsv.parquet.gz")

setDT(df_cells)
setDT(df_mcp)

df_cells[df_cell_map, cell_label := factor(i.cell_type, df_cell_map$cell_type),
    on = c(cell_type = "wg2_scpred_prediction")
]
df_stats[df_cell_map, cell_label := factor(i.cell_type, df_cell_map$cell_type),
    on = c(cell_type = "wg2_scpred_prediction")
]

# Visualize overall cell stats / count
df_cell_crohn <- inner_join(
    df_cells %>% select(index, cell_type, norm_score = crohns),
    df_mcp %>% select(index, cell_type, mc_pval = crohns)
)

df_cell_crohn2 <- read_parquet("results/misc/tenk10k_prop_sample.crohns.cell_score.tsv.parquet.gz")

df_crohns <- bind_rows(
    list(
        `10K-max per cell type` = df_cell_crohn,
        `10% per cell type` = df_cell_crohn2
    ),
    .id = "sampling"
)
df_crohns[df_cell_map, cell_label := factor(i.cell_type, df_cell_map$cell_type),
    on = c(cell_type = "wg2_scpred_prediction")
]

# Visualize distribution of scores
(p_violin <- ggplot(df_crohns) +
    theme_bw() +
    geom_violin(
        aes(
            x = norm_score, y = fct_rev(cell_label),
            color = cell_label,
            fill = after_scale(alpha(colour, 0.5))
        ),
        linewidth = 0.5
    ) +
    geom_vline(aes(xintercept = 0.0), linetype = "dashed", colour = "red") +
    facet_grid(cols = vars(sampling), scale = "free", switch = "y") +
    labs(x = "scDRS Crohn's disease", y = NULL) +
    scale_color_manual(
        values = deframe(df_cell_map[, .(cell_type, color)]),
        name = NULL,
        guide = guide_legend(ncol = 1)
    ) +
    theme(
        strip.background = element_blank(),
        legend.key.spacing.y = unit(0.25, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.5, "lines"),
        strip.text.x = element_text(face = "bold")
    )
)


# Plot proportion of enriched cells
df_crohns_sum <- df_crohns %>%
    mutate(`MC P-value < 0.05` = mc_pval < 0.05) %>%
    group_by(sampling, cell_label, `MC P-value < 0.05`) %>%
    tally() %>%
    mutate(prop = n / sum(n))

(p_compare_prop <- ggplot(df_crohns_sum) +
    theme_bw() +
    geom_col(
        aes(
            x = prop, y = fct_rev(cell_label),
            fill = `MC P-value < 0.05`, group = `MC P-value < 0.05`
        ),
        colour = "black", linewidth = 0.5, width = 1
    ) +
    facet_grid(cols = vars(sampling), scale = "free", switch = "y") +
    labs(x = "Cell proportion", y = NULL) +
    scale_fill_manual(
        values = c(`TRUE` = "red3", `FALSE` = "gray90"),
        breaks = c("TRUE", "FALSE"),
        labels = list(
            bquote(italic(P)[scDRS] < 0.05),
            bquote(italic(P)[scDRS] >= 0.05)
        ),
        name = NULL,
        guide = guide_legend(ncol = 1)
    ) +
    scale_x_continuous(labels = scales::percent) +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        strip.clip = "off",
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        legend.key.spacing.y = unit(0.25, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.5, "lines"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    )
)

ggsave("figures/supp/10-scdrs_crohns_sampling_violin.png",
    p_violin,
    width = 10, height = 8, device = agg_png
)

ggsave("figures/supp/11-scdrs_crohns_sampling_prop.png",
    p_compare_prop,
    width = 10, height = 8, device = agg_png
)
