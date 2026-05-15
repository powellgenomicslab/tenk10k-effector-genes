# smr overview plot
source("scripts/0-preprocess/preprocess_results.R")

library(scales)
library(patchwork)
library(ragg)
library(paletteer)
library(geomtextpath)

# cell type summary
df_cell_summary <- df_msmr_tenk10k[,
    .(
        n_genes = length(unique(probeID)),
        n_mr_genes = .SD[sig == TRUE, length(unique(probeID))],
        n_sig = sum(sig),
        n_tested = .N,
        prop_sig = sum(sig) / .N
    ),
    by = list(cell_type, major_cell_type)
]

(p_scatter_cell <- df_cell_summary %>%
    ggplot(aes(x = n_sig, y = n_tested)) +
    theme_bw() +
    geom_smooth(color = "gray20", linewidth = 0.5) +
    geom_point(aes(color = cell_type)) +
    scale_color_manual(
        values = deframe(df_cell_map[, .(cell_type, color)]),
        name = NULL
    ) +
    labs(
        x = "Number of significant MR associations",
        y = "Number of MR tests"
    ) +
    scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
    theme(
        aspect.ratio = 1,
        legend.text = element_text(size = 9)
    )
)

ggsave("figures/supp/3-smr_celltype_scatter.png",
    p_scatter_cell,
    scaling = 1.5,
    width = 10, height = 6, device = agg_png
)

df_summary_gene_pheno <- df_msmr_tenk10k %>%
    group_by(pheno_label, pheno_cat, probeID) %>%
    summarise(mr = max(sig), gwas = max(magma_gene)) %>%
    summarise(
        n_mr = sum(mr),
        n_gwas = sum(gwas),
        prop_mr_gwas = sum(mr & gwas) / sum(gwas),
        prop_gwas_mr = sum(mr & gwas) / sum(mr),
        mr_gwas_ratio = sum(mr) / sum(gwas),
        intersect_mr_gwas = sum(mr & gwas),
        mr_only = sum(mr & !gwas),
        gwas_only = sum(gwas & !mr),
        union_mr_gwas = sum(mr | gwas),
        jaccard_mr_gwas = sum(mr & gwas) / sum(mr | gwas)
    ) %>%
    ungroup() %>%
    mutate(pheno_label = fct_reorder(pheno_label, n_mr))

setDT(df_summary_gene_pheno)

df_summary_gene_pheno[df_trait_map, n_eff_gwas := i.n_eff,
    on = c(pheno_label = "label")
]

# scatter plot
rho_gwas_mr <- cor.test(df_summary_gene_pheno$n_gwas, df_summary_gene_pheno$n_mr)

(p_cor_gwas_mr <- ggplot(
    df_summary_gene_pheno,
    aes(x = n_gwas, y = n_mr)
) +
    theme_bw() +
    geom_point(aes(color = pheno_cat), size = 3, alpha = 0.8) +
    geom_smooth(color = "gray20", linewidth = 0.5) +
    scale_color_manual(
        values = deframe(df_trait_cat[, .(cat_order, color)]),
        name = "Phenotype category"
    ) +
    labs(
        x = "Number of GWAS genes",
        y = "Number of MR genes"
    ) +
    scale_x_continuous(expand = 0.01) +
    scale_y_continuous(expand = 0.01) +
    geom_abline(slope = 1, intercept = 0, color = "gray50", linewidth = 0.5, linetype = "dashed") +
    coord_equal() +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        # aspect.ratio = 1,
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0, 1)
    )
)

ggsave(p_cor_gwas_mr,
    filename = "figures/supp/4-smr_gwas_correlation.png",
    width = 6, height = 8, device = agg_png, scaling = 1.2
)

# get numbers
df_summary_gene_pheno[, list(
    n_trait = .N,
    n_mr_trait = sum(n_mr),
    avg_mr = mean(n_mr)
), ,
by = pheno_cat
]

# Comparison between MR genes - GWAS genes and TenK10K MR genes - eQTLgen

df_gene_summary_by_pheno <- df_msmr %>%
    rename(gene = probeID) %>%
    group_by(gene, phenotype, pheno_label, pheno_cat) %>%
    summarise(
        n_celltypes = n(),
        gwas = max(magma_gene, na.rm = TRUE) %>% as.logical(),
        eqtlgen = max(eqtlgen_mr, na.rm = TRUE) %>% as.logical()
    ) %>%
    pivot_longer(c(gwas, eqtlgen),
        names_to = "annotation", values_to = "value"
    ) %>%
    group_by(phenotype, pheno_label, pheno_cat, annotation, value) %>%
    tally(name = "n_genes") %>%
    ungroup() %>%
    complete(nesting(phenotype, pheno_label, pheno_cat, annotation),
        nesting(value),
        fill = list(n_genes = 0)
    ) %>%
    group_by(phenotype, pheno_label, pheno_cat, annotation) %>%
    mutate(
        n_mr_genes = sum(n_genes),
        prop_genes = n_genes / sum(n_genes)
    ) %>%
    mutate(fct_pheno = fct_reorder(pheno_label, n_mr_genes, mean) %>% fct_rev())

pals <- c("TRUE" = "#4C6C94FF", "FALSE" = "#A4ABB0FF")

mkplot_gene_summary <- function(df, x_col, y_col = fct_pheno,
                                xlab = NULL, col_scheme = pals,
                                labs = waiver(), name = NULL) {
    df %>%
        ggplot(aes(y = {{ y_col }}, x = {{ x_col }})) +
        theme_bw() +
        facet_wrap(~pheno_cat, nrow = 1, scale = "free_y", space = "free_y") +
        geom_col(aes(fill = value), width = 1, color = "black", linewidth = 0.4) +
        scale_fill_manual(
            values = pals, breaks = names(pals),
            labels = labs, name = NULL
        ) +
        labs(y = NULL, x = xlab) +
        scale_x_continuous(expand = expansion(0)) +
        theme(
            strip.text = element_text(
                hjust = 0, face = "bold", size = 9,
                margin = margin(b = 0.2, unit = "lines")
            ),
            strip.background = element_blank(),
            legend.key.size = unit(0.5, "lines"),
            legend.text = element_text(size = 9, margin = margin(l = 2)),
            legend.key.spacing.y = unit(0.25, "lines"),
            axis.title.x = element_text(size = 10),
            axis.text.x = element_text(size = 9),
            # panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.2, "lines")
        )
}


(p_gwas_n <- df_gene_summary_by_pheno %>%
    filter(annotation == "gwas") %>%
    mkplot_gene_summary(
        x_col = n_genes,
        xlab = "Number of unique MR genes",
        labs = c(
            "FALSE" = "MR only",
            "TRUE" = "MR & GWAS"
        )
    ) +
    theme(
        plot.margin = margin(),
        legend.position = "none"
    )
)

(p_gwas_prop <- df_gene_summary_by_pheno %>%
    filter(annotation == "gwas") %>%
    mkplot_gene_summary(
        x_col = prop_genes,
        xlab = "Proportion of MR genes",
        labs = c(
            "FALSE" = "MR only",
            "TRUE" = "MR & GWAS"
        )
    ) +
    scale_x_continuous(labels = percent, expand = expansion(0)) +
    geom_label(
        aes(
            label = ifelse(value == TRUE,
                paste0(" ", number(n_genes, 1, big.mark = ",")),
                paste0(number(n_genes, 1, big.mark = ","), " ")
            ),
            x = ifelse(value == TRUE, -Inf, Inf),
            hjust = ifelse(value == TRUE, 0, 1)
        ),
        linewidth = 0, label.padding = unit(0, "lines"),
        fill = alpha("white", 0.5),
        vjust = 0.5, size = 7 / .pt,
        position = position_nudge(x = 0.05)
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(color = NA),
        axis.ticks.y = element_line(color = "gray"),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 0),
        plot.margin = margin(r = 0.5, unit = "lines"),
        axis.ticks.length.y = unit(1, "lines")
    )
)

p_gwas_mr <- p_gwas_n + p_gwas_prop +
    plot_layout(nrow = 1, guides = "keep") &
    layer_scales(p_gwas_n)$y

ggsave("figures/supp/5-gwas_mr_gene_summary.png",
    p_gwas_mr,
    scaling = 1.05,
    width = 8, height = 12, device = agg_png
)

# eQTLgen MR genes
(p_eqtlgen_n <- df_gene_summary_by_pheno %>%
    filter(annotation == "eqtlgen") %>%
    mkplot_gene_summary(
        x_col = n_genes,
        xlab = "Number of unique MR genes",
        labs = c(
            "FALSE" = "TenK10K MR only",
            "TRUE" = "TenK10K & eQTLgen MR"
        )
    ) +
    theme(
        plot.margin = margin(),
        legend.position = "none"
    )
)

(p_eqtlgen_prop <- df_gene_summary_by_pheno %>%
    filter(annotation == "eqtlgen") %>%
    mkplot_gene_summary(
        x_col = prop_genes,
        xlab = "Proportion of MR genes",
        labs = c(
            "FALSE" = "TenK10K MR only",
            "TRUE" = "TenK10K & eQTLgen MR"
        )
    ) +
    scale_x_continuous(labels = percent, expand = expansion(0)) +
    geom_label(
        aes(
            label = ifelse(value == TRUE,
                paste0(" ", number(n_genes, 1, big.mark = ",")),
                paste0(number(n_genes, 1, big.mark = ","), " ")
            ),
            x = ifelse(value == TRUE, -Inf, Inf),
            hjust = ifelse(value == TRUE, 0, 1)
        ),
        linewidth = 0, label.padding = unit(0, "lines"),
        fill = alpha("white", 0.5),
        vjust = 0.5, size = 7 / .pt,
        position = position_nudge(x = 0.05)
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(color = NA),
        axis.ticks.y = element_line(color = "gray"),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 0),
        plot.margin = margin(r = 0.5, unit = "lines"),
        axis.ticks.length.y = unit(1, "lines")
    )
)


(p_eqtlgen <- p_eqtlgen_n + p_eqtlgen_prop +
    plot_layout(nrow = 1, guides = "keep") &
    layer_scales(p_eqtlgen_n)$y
)

ggsave("figures/supp/6-eqtlgen_mr_gene_summary.png",
    p_eqtlgen,
    scaling = 1.05,
    width = 8, height = 12, device = agg_png
)

# get number:
setDT(df_gene_summary_by_pheno)
tab_gene_summary_by_pheno <- df_gene_summary_by_pheno %>%
    pivot_wider(
        names_from = c(annotation, value), values_from = c(n_genes, prop_genes),
        names_sep = "."
    ) %>%
    setDT()

# write to supplementary table
source("scripts/util/write_table.R")
write_table(tab_gene_summary_by_pheno, "mr_gwas_eqtlgen", 4)
