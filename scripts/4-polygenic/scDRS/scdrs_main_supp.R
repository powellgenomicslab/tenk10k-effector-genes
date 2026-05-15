# Scripts to produce scDRS display items (main figure, supplementary figures, supplementarytables)

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
library(glue)
library(geomtextpath)

df_cell_map <- fread("metadata/cell.tsv") %>%
    mutate(cell_type = factor(cell_type, unique(cell_type)))

df_trait_map <- fread("metadata/trait.tsv")[include == TRUE]

df_trait_cat <- fread("metadata/trait_category.tsv")

df_stats <- fread("results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv") %>%
    filter(phenotype %in% df_trait_map$trait_id)
df_top <- fread("results/aggregate/tenk10k_phase1.scdrs.cell_type_top.tsv")

df_cell_map <- fread("resources/misc/cell_map.tsv") %>%
    mutate(cell_type = factor(cell_type, unique(cell_type)))
df_cells <- read_parquet("results/aggregate/tenk10k_phase1.scdrs.cell_score.tsv.parquet.gz")

# get cell mcp value
df_mcp <- read_parquet("results/aggregate/tenk10k_phase1.scdrs.cell_mcp.tsv.parquet.gz") %>%
    pivot_longer(any_of(df_trait_map$trait_id), names_to = "phenotype", values_to = "assoc_mcp") %>%
    setDT()

n_cells <- nrow(df_cells)


setDT(df_cells)
df_cells[df_cell_map,
    `:=`(
        cell_label = factor(i.cell_type, df_cell_map$cell_type),
        major_cell_type = factor(i.major_cell_type, unique(df_cell_map$major_cell_type))
    ),
    on = c(cell_type = "wg2_scpred_prediction")
]
df_stats[df_cell_map,
    `:=`(
        cell_label = factor(i.cell_type, df_cell_map$cell_type),
        major_cell_type = factor(i.major_cell_type, unique(df_cell_map$major_cell_type))
    ),
    on = c(cell_type = "wg2_scpred_prediction")
]
df_mcp[df_cell_map,
    `:=`(
        cell_label = factor(i.cell_type, df_cell_map$cell_type),
        major_cell_type = factor(i.major_cell_type, unique(df_cell_map$major_cell_type))
    ),
    on = c(cell_type = "wg2_scpred_prediction")
]

df_stats[df_trait_map, `:=`(
    pheno_label = i.label,
    pheno_cat = i.cat_rev,
    pheno_slim = i.include,
    supercategory = i.supercategory
),
on = c(phenotype = "trait_id")
]

df_mcp[df_trait_map, `:=`(
    pheno_label = i.label,
    pheno_cat = i.cat_rev,
    pheno_slim = i.include,
    supercategory = i.supercategory
),
on = c(phenotype = "trait_id")
]
df_stats[, pheno_cat := factor(pheno_cat, df_trait_cat$cat_order)]
df_mcp[, pheno_cat := factor(pheno_cat, df_trait_cat$cat_order)]

df_mcp_sig <- df_mcp[assoc_mcp < 0.05]

df_mcp_sig[df_stats, `:=`(
    ct_sig = i.assoc_mcp < 0.05,
    n_cell_ct = i.n_cell
),
on = c("phenotype", "cell_label")
]

df_mcp_sig_count <- df_mcp_sig %>%
    group_by(cell_label, major_cell_type, ct_sig, n_cell_ct, pheno_label, pheno_cat) %>%
    tally() %>%
    group_by(pheno_label, pheno_cat) %>%
    mutate(
        cell_rank = frank(n, ties.method = "first"),
        prop = n * as.numeric(ct_sig) / sum(as.numeric(ct_sig) * n_cell_ct)
    ) %>%
    arrange(pheno_cat, pheno_label, -prop) %>%
    mutate(cum_prop = cumsum(prop), lag_cum_prop = lag(cum_prop, default = 0)) %>%
    setDT()

df_mcp_sig_count[df_stats, `:=`(
    assoc_mcp = i.assoc_mcp,
    n_cell = i.n_cell,
    assoc_mcz = i.assoc_mcz
),
on = c("pheno_label", "cell_label")
]

df_pheno_order <- df_mcp_sig_count %>%
    group_by(pheno_label, pheno_cat) %>%
    summarise(
        n = sum(n * as.numeric(assoc_mcp < 0.05) * as.numeric(ct_sig)),
        n_prop = n / sum(as.numeric(ct_sig) * n_cell_ct)
    ) %>%
    arrange(pheno_cat, desc(n_prop))

(p_sig_count <- df_mcp_sig_count %>%
    group_by(pheno_cat) %>%
    mutate(pheno_label = factor(pheno_label, df_pheno_order$pheno_label)) %>%
    ggplot(aes(x = prop, y = pheno_label, group = pheno_label)) +
    theme_void() +
    geom_col(aes(fill = cell_label),
        data = ~ filter(.x, assoc_mcp < 0.05),
        color = "black", linewidth = 0.5
    ) +
    facet_wrap(vars(pheno_cat),
        dir = "h", ncol = 1,
        space = "free_y", scales = "free_y"
    ) +
    geom_text(aes(label = paste0(percent(x, 1), " "), x = 0),
        hjust = 1,
        size = 10 / .pt,
        data = ~ .x %>%
            filter(ct_sig) %>%
            group_by(pheno_cat, pheno_label) %>%
            summarise(x = max(cum_prop))
    ) +
    geom_text(aes(label = pheno_label, x = x),
        hjust = 0,
        size = 11 / .pt, nudge_x = 0.01,
        data = ~ .x %>%
            filter(ct_sig) %>%
            group_by(pheno_cat, pheno_label) %>%
            summarise(x = max(cum_prop))
    ) +
    scale_fill_manual(
        values = deframe(df_cell_map[, .(cell_type, color)]),
        breaks = df_cell_map$cell_type,
        name = NULL,
        guide = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)
    ) +
    scale_x_continuous(
        expand = expansion(mult = c(0.01, 0.01)),
        labels = scales::percent,
        guide = guide_axis(cap = "lower")
    ) +
    scale_y_discrete(expand = expansion(add = c(1, 0.6))) +
    coord_cartesian(clip = "off") +
    labs(x = "Average proportion of phenotype-enriched cells") +
    theme(
        legend.position = "inside",
        panel.spacing.y = unit(0.1, "lines"),
        axis.text.x = element_text(size = 10, hjust = 0.5),
        axis.line.x = element_line(arrow = arrow(length = unit(0.2, "lines"), ends = "last")),
        axis.ticks.x = element_line(),
        axis.title.x = element_text(size = 10),
        axis.ticks.length.x = unit(0.25, "lines"),
        plot.margin = margin(r = 10, l = 1.5, unit = "lines"),
        legend.text = element_text(size = 9, margin = margin(l = 2)),
        legend.key.size = unit(0.5, "lines"),
        legend.key.spacing.y = unit(0.1, "lines"),
        legend.position.inside = c(1.2, 0),
        legend.justification = c("left", "bottom"),
        strip.text = element_text(
            hjust = 0, face = "bold", size = 12,
            margin = margin(b = 0.25, unit = "lines")
        ),
        strip.background = element_blank()
    )
)

# export to supplementary tables
source("scripts/util/write_table.R")
tab_scdrs <- df_stats %>%
    left_join(df_mcp_sig_count) %>%
    mutate(prop_ct = n / n_cell_ct) %>%
    mutate(pheno_label = factor(pheno_label, df_pheno_order$pheno_label)) %>%
    arrange(pheno_cat, pheno_label, desc(prop))

write_table(tab_scdrs, "scdrs_results", 5)

# Visualize overall cell stats / count

# number of enriched traits across cell types
df_n_traits_by_cells <- df_mcp_sig_count %>%
    group_by(cell_label, pheno_cat) %>%
    mutate(
        pheno_length = n(),
        pheno_rank = dense_rank(n)
    )

trait_cat_col <- select(df_trait_cat, cat_order, color) %>% deframe()

df_stats %>%
    ggplot(aes(x = pheno_cat, y = -log10(assoc_mcp))) +
    geom_point(
        data = ~ filter(.x, assoc_mcp >= 0.05),
        position = position_jitter(width = 0.2, height = 0),
        colour = "gray80", shape = "circle filled",
        fill = alpha("gray80", 0.5)
    ) +
    geom_point(aes(fill = pheno_cat),
        alpha = 0.8,
        shape = "circle filled", colour = "black",
        position = position_jitter(width = 0.2, height = 0),
        data = ~ filter(.x, assoc_mcp < 0.05)
    ) +
    theme_bw() +
    geom_hline(
        yintercept = -log10(0.05), linetype = "dashed",
        color = "red3"
    ) +
    labs(x = NULL, y = bquote(-log[10] ~ italic("P")["phenotype enrichment"])) +
    facet_wrap(~cell_label, scales = "free_x", nrow = 7) +
    scale_fill_manual(values = trait_cat_col, name = "Phenotype category") +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

# significance count of phenotypes by cell type
df_pheno_count <- df_mcp_sig_count %>%
    group_by(cell_label, major_cell_type, pheno_label, pheno_cat) %>%
    summarise(
        sig_count = sum(n * as.numeric(assoc_mcp < 0.05)),
        sig_prop = sum(n / n_cell * as.numeric(assoc_mcp < 0.05))
    ) %>%
    ungroup() %>%
    complete(nesting(pheno_label, pheno_cat), nesting(cell_label, major_cell_type),
        fill = list(sig_count = 0, sig_prop = 0)
    )

(p_grid <- df_pheno_count %>%
    ggplot(aes(y = pheno_label, x = cell_label)) +
    geom_tile(
        data = ~ filter(.x, sig_prop == 0),
        fill = "gray90", color = "white", linewidth = 0.2
    ) +
    geom_tile(aes(fill = sig_prop),
        data = ~ filter(.x, sig_prop > 0),
        color = "black", linewidth = 0.5
    ) +
    theme_bw() +
    scale_x_discrete(expand = expansion(add = c(0))) +
    scale_y_discrete(expand = expansion(add = c(0))) +
    facet_grid(
        scale = "free", space = "free",
        rows = vars(pheno_cat),
        cols = vars(major_cell_type)
    ) +
    scale_fill_paletteer_c("grDevices::Plasma",
        direction = -1,
        limits = c(0, 1),
        labels = scales::percent
    ) +
    labs(
        fill = "Proportion of phenotype-enriched cells",
        y = NULL, x = NULL
    ) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, margin = margin(b = 1, r = 1, unit = "lines")),
        legend.title.position = "left",
        plot.margin = margin(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)

(p_n_celltype_by_pheno <- df_pheno_count %>%
    group_by(pheno_label, pheno_cat) %>%
    summarise(sig_count = sum(sig_count > 0)) %>%
    ungroup() %>%
    mutate(mean = mean(sig_count)) %>%
    ggplot(aes(x = sig_count, y = pheno_label)) +
    facet_grid(
        scale = "free", space = "free",
        rows = vars(pheno_cat)
    ) +
    theme_void() +
    geom_col(
        fill = "steelblue3", color = "black", width = 1,
        linewidth = 0.5
    ) +
    scale_x_continuous(expand = expansion(add = c(0.0, 2))) +
    layer_scales(p_grid)$y +
    geom_text(aes(label = sig_count), nudge_x = 0.2, hjust = 0, size = 8 / .pt) +
    geom_text(aes(label = label, x = 5, y = -Inf),
        hjust = 0.5,
        size = 10 / .pt, vjust = 2,
        data = ~ .x %>%
            distinct(pheno_cat) %>%
            mutate(label = ifelse(pheno_cat == "Other", "N cell types", NA))
    ) +
    geom_vline(aes(xintercept = mean),
        linetype = "dashed",
        data = ~ distinct(.x, pheno_cat, mean),
        color = "red3", linewidth = 0.5
    ) +
    geom_label(aes(label = label, y = -Inf, x = mean),
        hjust = 0, vjust = 0,
        fill = "white", size = 10 / .pt, linewidth = 0,
        data = ~ distinct(.x, pheno_cat, mean) %>%
            mutate(label = ifelse(pheno_cat == "Other",
                paste("Mean:", number(mean, 0.1)), NA
            )),
        color = "red3", linewidth = 0.5
    ) +
    coord_cartesian(clip = "off") +
    theme(
        strip.text.y = element_text(
            face = "bold", hjust = 0,
            margin = margin(l = 0.2, unit = "lines")
        ),
        panel.border = element_rect(color = "black", linewidth = 0.5)
    )
)

(p_pheno_by_celltype <- df_pheno_count %>%
    group_by(cell_label, major_cell_type) %>%
    summarise(sig_count = sum(sig_count > 0)) %>%
    ungroup() %>%
    mutate(mean = mean(sig_count)) %>%
    ggplot(aes(y = sig_count, x = cell_label)) +
    facet_grid(
        scale = "free", space = "free",
        cols = vars(major_cell_type)
    ) +
    layer_scales(p_grid)$x +
    scale_y_continuous(expand = expansion(add = c(0.0, 4))) +
    geom_text(aes(label = sig_count), nudge_y = 0.5, size = 8 / .pt, vjust = 0) +
    geom_text(aes(label = label, x = -Inf, y = 20),
        hjust = 0.5,
        size = 10 / .pt, vjust = -1, angle = 90,
        data = ~ .x %>%
            distinct(major_cell_type) %>%
            mutate(label = ifelse(major_cell_type == "CD4 T",
                "N phenotypes", NA
            ))
    ) +
    theme_void() +
    coord_cartesian(clip = "off") +
    geom_col(fill = "steelblue3", color = "black", width = 1, linewidth = 0.5) +
    geom_hline(aes(yintercept = mean),
        linetype = "dashed",
        data = ~ distinct(.x, major_cell_type, mean),
        color = "red3", linewidth = 0.5
    ) +
    geom_label(aes(label = label, x = -Inf, y = mean),
        hjust = 0, vjust = 0,
        fill = "white", size = 10 / .pt, linewidth = 0,
        data = ~ distinct(.x, major_cell_type, mean) %>%
            mutate(label = ifelse(major_cell_type == "CD4 T",
                paste("Mean:", number(mean, 0.1)), NA
            )),
        color = "red3", linewidth = 0.5
    ) +
    theme(
        strip.text.x = element_text(
            face = "bold",
            margin = margin(b = 0.2, unit = "lines")
        ),
        strip.clip = "off",
        panel.border = element_rect(color = "black", linewidth = 0.5)
    )
)

p_blank <- ggplot() +
    theme_void() +
    layer_scales(p_n_celltype_by_pheno)$x +
    layer_scales(p_pheno_by_celltype)$y

(plots_grid <- list(
    p_pheno_by_celltype, p_blank,
    p_grid, p_n_celltype_by_pheno
) %>%
    map(~ .x + theme(panel.spacing = unit(0.5, "lines"))) %>%
    wrap_plots(
        heights = c(0.1, 1),
        widths = c(1, 0.2),
        guides = "keep"
    )
)


ggsave("figures/supp/7-scdrs_pheno_celltype_grid.png",
    plots_grid,
    bg = "white", width = 12, height = 16, device = agg_png, scaling = 1
)


# Crohn's disease example

# Visualize cell score UMAP for crohns
df_crohns <- df_cells %>%
    select(index, cell_type, cell_label, major_cell_type, umap_1, umap_2, scdrs = crohns) %>%
    mutate(
        across(c(umap_1, umap_2), ~ .x %in% boxplot.stats(.x)$out,
            .names = "outlier_{col}"
        ),
        .by = cell_type
    )

df_crohns[df_mcp[phenotype == "crohns"], assoc_mcp := i.assoc_mcp, on = "index"]
df_crohns[df_stats[phenotype == "crohns" & assoc_mcp < 0.05], ct_sig := TRUE, on = c("cell_type", "cell_label")]

# calculate centroid
centroid <- function(x) mean(x[!x %in% boxplot.stats(x)$out], na.rm = T)

df_labels <- df_crohns[ct_sig == TRUE, .(
    umap_1 = centroid(umap_1),
    umap_2 = centroid(umap_2)
),
by = cell_label
]

set.seed(100)
(p_scdrs_umap_crohns <- df_crohns %>%
    slice_sample(n = nrow(.)) %>%
    # arrange(-assoc_mcp) %>%
    ggplot(aes(x = umap_1, y = umap_2)) +
    theme_void() +
    # geom_point(data = ~filter(.x, is.na(assoc_mcp) | assoc_mcp >= 0.05),
    #            shape = "circle", colour = "#ECEFF1FF",
    #            size = 0.5, alpha = 0.5) +
    # geom_point(aes(colour = ifelse(is.na(ct_sig), NA_character_, as.character(cell_label))),
    geom_point(aes(colour = ifelse(is.na(ct_sig), NA, scdrs)),
        # data = ~filter(.x, assoc_mcp < 0.05 & ct_sig),
        shape = "circle",
        # show.legend = TRUE,
        size = 0.2, alpha = 0.5
    ) +
    geom_label_repel(
        aes(label = cell_label),
        color = "black",
        fontface = "bold", force = 1, label.size = NA, label.padding = 0.1,
        size = 12 / .pt, fill = alpha("white", 0.8),
        data = df_labels, max.overlaps = Inf
    ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed() +
    scale_colour_viridis_c(
        name = "scDRS\nCrohn's\ndisease", option = "magma", direction = -1,
        na.value = "gray90"
    ) +
    theme(
        panel.border = element_rect(linewidth = 1, color = "black"),
        legend.position = "right",
        legend.position.inside = c(0.75, 0.02),
        axis.title.x = element_text(size = 12, margin = margin(t = 0.5, unit = "lines")),
        axis.title.y = element_text(
            size = 12, margin = margin(r = 0.5, unit = "lines"),
            angle = 90
        )
    )
)

# Ridge plot
df_cell_order <- df_stats[phenotype == "crohns"] %>%
    arrange(-assoc_mcz) %>%
    mutate(
        cell_label = factor(cell_label, unique(cell_label)),
        sig = assoc_mcp < 0.05,
        text_col = ifelse(sig, "red3", "black"),
        text_face = ifelse(sig, "bold", "plain")
    )
df_avg <- df_crohns[, .(
    mean_scdrs = mean(scdrs, na.rm = TRUE),
    sd_scdrs = sd(scdrs, na.rm = TRUE),
    n_cells = .N,
    ct_sig = !is.na(max(ct_sig)),
    n_sig = sum(as.numeric(assoc_mcp < 0.05)),
    prop_sig = sum(as.numeric(assoc_mcp < 0.05)) / .N
),
by = cell_label
]

plusmin <- "±"
(p_ridge <- df_crohns %>%
    ggplot(aes(x = scdrs, y = factor(cell_label, rev(df_cell_order$cell_label)))) +
    theme_bw() +
    geom_density_ridges(
        aes(colour = !is.na(ct_sig), fill = after_scale(alpha(colour, 0.5))),
        scale = 1.4, linewidth = 0.2, rel_min_height = 0,
    ) +
    labs(x = "scDRS Crohn's disease", y = NULL) +
    scale_colour_manual(
        values = c("FALSE" = "gray50", "TRUE" = "red3"),
        name = NULL,
        guide = "none"
    ) +
    geom_text(aes(label = paste(number(mean_scdrs, 0.1), plusmin, number(sd_scdrs, 0.1))),
        data = df_avg[ct_sig == TRUE], size = 10 / .pt, color = "red3",
        x = Inf, vjust = 0, position = position_nudge(y = 0.1), hjust = 1
    ) +
    geom_text(aes(label = paste(number(mean_scdrs, 0.1), plusmin, number(sd_scdrs, 0.1))),
        data = df_avg[ct_sig == FALSE], size = 10 / .pt, color = "black",
        x = Inf, vjust = 0, position = position_nudge(y = 0.1), hjust = 1
    ) +
    annotate("text",
        x = Inf, y = Inf, hjust = 1, vjust = 1,
        fontface = "plain", label = paste("Mean", plusmin, "SD"), size = 10 / .pt
    ) +
    scale_x_continuous(
        labels = label_number(accuracy = 0.01),
        expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme(
        axis.title.x = element_text(size = 11),
        axis.text.y = element_text(
            vjust = 0, size = 11,
            face = rev(df_cell_order$text_face),
            colour = rev(df_cell_order$text_col)
        ),
        axis.text.x = element_text(size = 11),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.3, "lines"),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        plot.margin = margin(l = -5, unit = "lines"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank()
    )
)

(p_bar_prop <- df_avg %>%
    ggplot(aes(y = factor(cell_label, rev(df_cell_order$cell_label)))) +
    geom_col(aes(x = n_cells),
        position = position_nudge(y = 0.25),
        fill = "gray80", width = 0.5,
        color = "black", linewidth = 0.5
    ) +
    geom_col(aes(x = n_sig, fill = "red3"),
        position = position_nudge(y = 0.25),
        # fill = "#4E4E4EFF",
        colour = "black", width = 0.5, linewidth = 0.5
    ) +
    theme_bw() +
    scale_fill_identity(
        name = NULL, labels = expression(italic("P")["scDRS"] < 0.05),
        guide = guide_legend()
    ) +
    geom_text(aes(x = Inf, color = ct_sig, label = paste0(" ", percent(prop_sig, 0.1))),
        size = 10 / .pt, position = position_nudge(y = 0.1),
        hjust = 1, vjust = 0
    ) +
    scale_color_manual(
        values = c("FALSE" = "black", "TRUE" = "red3"),
        guide = "none"
    ) +
    labs(x = "Number of cells", y = NULL) +
    annotate("text",
        x = Inf, y = Inf, hjust = 1, vjust = 1,
        fontface = "plain", label = bquote("% enriched cells"), size = 10 / .pt
    ) +
    scale_x_continuous(
        expand = expansion(mult = c(0.01)),
        limits = c(0, 1.7e4),
        breaks = c(0, 5000, 10000),
        labels = c("0", "5k", "10k"),
        guide = guide_axis(cap = "both")
    ) +
    layer_scales(p_ridge)$y +
    coord_cartesian(clip = "off") +
    theme(
        axis.title.x = element_text(size = 11, hjust = 0),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.3, "lines"),
        axis.line.x = element_line(),
        panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(l = 0.5, unit = "lines"),
        legend.position = "none",
        panel.grid.major.y = element_blank()
    )
)


(p_ridge_bar <- (p_ridge + p_bar_prop) +
    plot_layout(width = c(1, 0.4))
)

# Decile plot
df_q_scdrs <- df_cells %>%
    select(crohns, B2_ALL_eur_leave_23andme, cell_label) %>%
    pivot_longer(-cell_label, names_to = "phenotype", values_to = "scdrs") %>%
    left_join(df_trait_map %>% select(phenotype = trait_id, pheno_label = name)) %>%
    group_by(pheno_label, phenotype, cell_label) %>%
    mutate(q_scdrs = ntile(scdrs, 10)) %>%
    group_by(pheno_label, phenotype, cell_label, q_scdrs) %>%
    summarise(mean = mean(scdrs), sd = sd(scdrs), n = n()) %>%
    setDT()

df_q_scdrs[df_stats, ct_sig := i.assoc_mcp < 0.05, on = c("cell_label", "phenotype")]

dendritic_cells <- c("ASDC", "cDC1", "cDC2", "pDC")

(p_decile <- df_q_scdrs %>%
    filter(cell_label %in% dendritic_cells) %>%
    ggplot(aes(
        x = q_scdrs, y = mean,
        alpha = ifelse(ct_sig, 1, 0.4),
        colour = pheno_label, group = pheno_label
    )) +
    facet_wrap(~cell_label, axes = "all") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 2) +
    geom_line(aes(label = pheno_label), linewidth = 1) +
    scale_alpha_identity() +
    labs(
        y = "Mean scDRS",
        x = "scDRS decile"
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(
        breaks = 1:10, labels = 1:10,
        limits = c(1, 10),
        expand = expansion(add = c(0.4)),
        guide = guide_axis(cap = "both")
    ) +
    scale_colour_paletteer_d(
        "ggthemes::Tableau_10",
        name = NULL
    ) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.2, "lines"),
        plot.margin = margin(l = 1, unit = "lines"),
        axis.line = element_line(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.ticks.length = unit(0.3, "lines")
    )
)

# combine plots

design <- "
AAAAABBBBB
AAAAABBBBB
AAAAACCCCC
AAAAACCCCC
AAAAACCCCC
AAAAADDDDD
AAAAADDDDD
"
list_p <- list(
    A = p_sig_count,
    B = p_scdrs_umap_crohns,
    C = p_ridge_bar,
    D = p_decile
) %>%
    map(wrap_elements)

(plots <- wrap_plots(list_p, design = design) +
    plot_annotation(tag_levels = "a") &
    theme(
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.01, 1)
    )
)

ggsave("figures/main/3-polygenic_enrichment.pdf",
    plots,
    bg = "white", width = 12.5, height = 18,
    device = cairo_pdf, scale = 1 / 1.2
)
