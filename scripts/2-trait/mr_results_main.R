# Fig 2: MR overview
source("scripts/0-preprocess/preprocess_results.R")

library(patchwork)
library(ragg)
library(scales)
library(paletteer)
library(ggnewscale)
library(geomtextpath)


# panel a1: tally by N celltypes and GWAS
# panel a2: tally by N celltypes and eQTLgen MR

mk_plot_tally <- function(group_col, label = list(mr_only = "MR only", mr_other = "MR & GWAS"),
                          p_title = NULL,
                          col_scheme = list(mr_only = "#A4ABB0FF", mr_other = "#4C6C94FF")) {
    df_tally <- df_msmr %>%
        group_by(probeID, phenotype, {{ group_col }}) %>%
        tally(name = "n_celltypes") %>%
        group_by(n_celltypes, {{ group_col }}) %>%
        tally(name = "n") %>%
        ungroup() %>%
        mutate(sum_n = sum(n), prop = n / sum_n, .by = n_celltypes) %>%
        mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))

    ggplot(df_tally, aes(x = n_celltypes, y = n)) +
        geom_col(aes(fill = {{ group_col }}),
            width = 1, linewidth = 0.5,
            position = "stack", color = "black"
        ) +
        scale_fill_manual(
            values = c(col_scheme$mr_only, col_scheme$mr_other),
            name = NULL,
            label = c(
                "FALSE" = label$mr_only,
                "TRUE" = label$mr_other
            ),
            guide = guide_legend(ncol = 1)
        ) +
        theme_classic() +
        labs(
            y = "N gene-trait MR associations",
            x = "N cell types with MR associations",
            title = p_title
        ) +
        geom_line(aes(group = {{ group_col }}, color = {{ group_col }}, y = prop_scaled),
            data = ~ filter(.x, !{{ group_col }})
        ) +
        geom_textline(
            aes(
                group = {{ group_col }},
                label = paste("%", ifelse({{ group_col }}, label$mr_other, label$mr_only)),
                color = {{ group_col }}, y = prop_scaled
            ),
            data = ~ filter(.x, !{{ group_col }}),
            size = 9 / .pt,
            vjust = -0.5, hjust = 0.2
        ) +
        geom_point(aes(color = {{ group_col }}, y = prop_scaled), data = ~ filter(.x, !{{ group_col }})) +
        scale_color_manual(
            values = c("red3", "red3"),
            guide = "none"
        ) +
        scale_x_continuous(
            breaks = seq(1, 28, 2),
            expand = expansion(add = 0.25)
        ) +
        scale_y_continuous(
            labels = ~ ifelse(.x < 1e3, .x, number(.x / 1e3, 1, suffix = "k")),
            n.breaks = 10,
            expand = expansion(mult = c(0, 0.05)),
            sec.axis = sec_axis(
                ~ rescale(., to = c(0, 1)),
                name = NULL,
                labels = label_percent()
            ),
            guide = guide_axis(cap = "none")
        ) +
        coord_cartesian(clip = "off") +
        theme(
            legend.position.inside = c(0.95, 1),
            legend.position = "inside",
            plot.title = element_text(size = 10, hjust = 0),
            axis.title = element_text(size = 9),
            legend.text = element_text(margin = margin(l = 2)),
            legend.justification = c(1, 1),
            legend.key.size = unit(0.75, "lines"),
            legend.key.spacing.y = unit(0.5, "lines"),
            axis.text = element_text(size = 8),
            # plot.margin = margin(t = 2, unit = "lines"),
            axis.line.y.right = element_line(color = "red3"),
            axis.text.y.right = element_text(color = "red3"),
            axis.ticks.y.right = element_line(color = "red3"),
            axis.title.y.right = element_text(color = "red3"),
            panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
        )
}

(p_mr_gwas <- mk_plot_tally(magma_gene, p_title = "TenK10K single-cell MR vs. GWAS"))
(p_mr_eqtlgen <- mk_plot_tally(eqtlgen_mr,
    p_title = "TenK10K single-cell MR vs. eQTLgen bulk whole blood MR",
    label = list(mr_only = "TenK10K only", mr_other = "TenK10K & eQTLgen")
)
)

# panel b: bar chart of direction of effect
two_color_pal <- paletteer::paletteer_d("RColorBrewer::RdBu")[c(3, 9)] %>%
    set_names(c("Positive", "Negative"))

df_msmr_by_cat_by_doe <- df_msmr %>%
    mutate(`Direction of effect` = fifelse(b_SMR > 0, "Positive", "Negative")) %>%
    mutate(`Cell type` = str_replace(biosample, "_", " ")) %>%
    count(`Cell type`, pheno_cat, supercategory, `Direction of effect`, name = "Number of SMR+HEIDI significant trait-gene pairs") %>%
    mutate(cat_totals = sum(`Number of SMR+HEIDI significant trait-gene pairs`), .by = c(pheno_cat, `Direction of effect`)) %>%
    mutate(`Cell type` = factor(`Cell type`, levels = df_cell_map$cell_type)) %>%
    rename("Phenotype category" = pheno_cat)

(p_bar_doe <- df_msmr_by_cat_by_doe %>%
    # ggplot(aes(x = `Cell type`, y = `Number of SMR+HEIDI significant trait-gene pairs`, fill = `Phenotype category`, alpha = `Direction of effect`)) +
    ggplot(aes(x = `Cell type`, y = `Number of SMR+HEIDI significant trait-gene pairs`, fill = `Direction of effect`)) +
    facet_wrap(~`Phenotype category`, ncol = 2, scales = "free_y") +
    geom_col() +
    labs(x = NULL, y = "N gene-trait MR associations") +
    scale_y_continuous(
        expand = expansion(0),
        labels = ~ ifelse(.x < 1e4, .x, label_number(1, scale = 1e-3, suffix = "k")(.x))
    ) +
    scale_fill_manual(
        values = two_color_pal,
        labels = c("Negative association", "Positive association"),
        name = NULL
    ) +
    # scale_alpha_manual(values = setNames(c(0.5, 1), c("Negative", "Positive"))) +
    guides(
        y = guide_axis(cap = "both"),
        fill = guide_legend(
            theme = theme(legend.key.size = unit(0.5, "lines")),
            direction = "horizontal"
        )
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, size = 7),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major.y = element_line(color = "gray90"),
        legend.position = "bottom",
        legend.text = element_text(margin = margin(l = 2)),
        # legend.position.inside = c(0.5, 1),
        # legend.justification.inside = c(1, 0),
        axis.title = element_text(size = 9),
        strip.text = element_text(size = 9, hjust = 0),
    )
)

# panel c: cell type by phenotype bar charts
trait_cat_col <- deframe(df_trait_cat[, .(cat_order, color)])

df_unique_gene <- df_msmr %>%
    group_by(pheno_label, pheno_cat) %>%
    summarise(n = n_distinct(probeID)) %>%
    ungroup() %>%
    filter(n >= 10) %>%
    mutate(fct_pheno = fct_reorder(pheno_label, -n)) %>%
    group_by(pheno_cat) %>%
    slice_max(n, n = 10) %>%
    setDT()

df_smr_pheno_cell <- df_msmr_tenk10k %>%
    filter(pheno_label %in% df_unique_gene$fct_pheno) %>%
    group_by(pheno_label, pheno_cat, supercategory, cell_type, major_cell_type) %>%
    summarise(n_mr = sum(sig, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(fct_pheno = factor(pheno_label, levels(df_unique_gene$fct_pheno))) %>%
    setDT()

(p_mr_summary <- df_smr_pheno_cell %>%
    ggplot(aes(y = n_mr, x = cell_type, fill = cell_type)) +
    facet_wrap(~fct_pheno, scale = "free_y", ncol = 6) +
    geom_col(color = "black", linewidth = 0.2, width = 1) +
    theme_bw() +
    scale_fill_manual(
        values = deframe(df_cell_map[, .(cell_type, color)]),
        name = "Cell type",
        guide = guide_legend(
            ncol = 1, theme = theme(
                legend.key.spacing.y = unit(0.15, "lines"),
                legend.title = element_text(size = 9, margin = margin(b = 4)),
                legend.key.size = unit(0.5, "lines")
            )
        )
    ) +
    scale_y_continuous(
        breaks = ~ unique(floor(pretty(seq(min(.x), max(.x) * 1.01), n = 4))),
        # labels = ~ifelse(.x < 1e3, .x, number(.x / 1e3, 0.1, suffix = "k")),
        expand = expansion(mult = c(0.02, 0.4))
    ) +
    labs(y = "Number of MR genes", x = NULL) +
    coord_cartesian(clip = "off") +
    ggnewscale::new_scale_fill() +
    geom_label(aes(label = n, color = "black", y = Inf, x = Inf),
        inherit.aes = FALSE,
        hjust = 1, vjust = 1,
        size = 7 / .pt, label.padding = unit(0.25, "lines"), label.size = 0,
        data = df_unique_gene
    ) +
    scale_color_identity(
        labels = "Number of\nunique\nMR genes",
        name = NULL,
        guide = guide_legend(
            order = 1,
            override.aes = list(label = "N", label.size = 0.5, vjust = 0.5, hjust = 0),
            theme = theme(legend.text = element_text(size = 8, vjust = 0.5, hjust = 0, margin = margin(l = -2)))
        )
    ) +
    geom_label(
        aes(
            label = n, fill = pheno_cat, y = Inf, x = Inf,
            color = after_scale(prismatic::best_contrast(fill, c("white", "black")))
        ),
        inherit.aes = FALSE, border.colour = "black",
        hjust = 1, vjust = 1,
        size = 7 / .pt, label.padding = unit(0.25, "lines"), label.size = 0.25,
        data = df_unique_gene
    ) +
    scale_fill_manual(
        values = trait_cat_col,
        name = "Trait category",
        guide = guide_legend(
            override.aes = list(label = "  ", vjust = 0.5, hjust = 0.5, size = 7 / .pt, lineheight = 0.5),
            theme = theme(
                legend.key.spacing.y = unit(-0.25, "lines"),
                legend.text = element_text(lineheight = 0.5, size = 8, vjust = 0.5, hjust = 0, margin = margin(l = -2))
            )
        )
    ) +
    theme(
        axis.text.x = element_blank(),
        strip.clip = "off",
        strip.background = element_blank(),
        strip.text = element_text(size = 8, margin = margin(t = 2, b = 2)),
        panel.spacing.y = unit(0.25, "lines"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.spacing.y = unit(0.25, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(2, "pt"),
        legend.title = element_text(size = 9, margin = margin(b = 1)),
        legend.box.spacing = unit(0.3, "lines"),
        legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0),
        # plot.margin = margin(t = 2, b = 2, l = 2, r = 0),
        # legend.key.spacing.y = unit(0, "lines"),
        legend.text = element_text(size = 8, vjust = 0.5, hjust = 0, margin = margin(l = 2, r = -5))
    )
)

# Combine

p_design <- "
AAABBBB
AAABBBB
DDDDDDD
DDDDDDD
DDDDDDD
"

(plots <- list(
    A = wrap_elements(full = p_mr_gwas / p_mr_eqtlgen),
    B = wrap_elements(full = p_bar_doe),
    D = wrap_elements(full = p_mr_summary)
) %>%
    wrap_plots(design = p_design, guides = "keep") +
    plot_annotation(
        tag_levels = "a",
        theme = theme(plot.margin = margin())
    ) &
    theme(
        plot.tag = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5),
        plot.tag.position = c(0.01, 0.99)
    )
)

ggsave("figures/main/2-mr_combined.pdf",
    dpi = 320,
    plots, width = 9, height = 12, device = cairo_pdf, scale = 1 / 0.9
)
