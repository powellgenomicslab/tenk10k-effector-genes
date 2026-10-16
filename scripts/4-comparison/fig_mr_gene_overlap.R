# figures: per-trait overlap between TenK10K MR genes and MAGMA / eQTLgen bulk-blood MR genes
#
# Recovered from scripts/mr_overview/overview_supplement.R (deleted in commit e137f97,
# 2026-05-12). Two figures:
#   - fig-s10 (gwas_mr_gene_summary):    TenK10K MR vs. MAGMA gene counts/proportions by trait
#   - fig-s11 (eqtlgen_mr_gene_summary): TenK10K MR vs. eQTLgen MR gene counts/proportions by trait
#
# The df_gene_summary_by_pheno build below matches scripts/table/tab_mr_overlap.R with
# one deliberate exception: this script omits the ungroup() before the fct_pheno mutate
# (see note at that line) so the trait ordering matches the originally published figure.

source("scripts/0-preprocess/preprocess_results.R")

library(patchwork)
library(ragg)
library(scales)
library(paletteer)

df_gene_summary_by_pheno <- df_msmr_tenk10k[
  probeID %in% gene_universe,
  .(mr = max(mr),
    magma_gene = max(magma_gene),
    eqtlgen_mr = max(eqtlgen_mr)),
  by = list(probeID, phenotype, pheno_label, pheno_cat)
] |>
  mutate(
    eqtlgen = case_when(
      eqtlgen_mr == 1 & mr == 1 ~ "TenK10K MR & eQTLgen MR",
      eqtlgen_mr == 0 & mr == 1 ~ "TenK10K MR only",
      eqtlgen_mr == 1 & mr == 0 ~ "eQTLgen MR only",
      TRUE ~ NA_character_),
    gwas = case_when(
      magma_gene == 1 & mr == 1 ~ "TenK10K MR & MAGMA",
      magma_gene == 0 & mr == 1 ~ "TenK10K MR only",
      magma_gene == 1 & mr == 0 ~ "MAGMA only",
      TRUE ~ NA_character_)
  ) |>
  pivot_longer(c(gwas, eqtlgen),
               values_drop_na = TRUE,
               names_to = "annotation", values_to = "value") |>
  group_by(phenotype, pheno_label, pheno_cat, annotation, value) |>
  tally(name = "n_genes") |>
  ungroup() |>
  complete(nesting(phenotype, pheno_label, pheno_cat),
           nesting(annotation, value),
           fill = list(n_genes = 0)) |>
  group_by(phenotype, pheno_label, pheno_cat, annotation) |>
  mutate(n_genes_total = sum(n_genes),
         prop_genes = n_genes / sum(n_genes)) |>
  # NOTE: no ungroup() here, matching the original deleted script exactly. fct_reorder()
  # below therefore runs inside the still-active group_by (each group is a single trait),
  # so it can't actually rank across traits -- the resulting order is an artifact of how
  # dplyr recombines the per-group single-level factors, not a deliberate sort. Kept as-is
  # to match the previously published figure's trait ordering.
  mutate(fct_pheno = fct_reorder(pheno_label, n_genes_total, mean) |> fct_rev()) |>
  setDT()

mkplot_gene_summary <- function(df, x_col, y_col = fct_pheno,
                                xlab = NULL, labs = waiver(), name = NULL) {
  df %>%
    ggplot(aes(y = {{y_col}}, x = {{x_col}})) +
    theme_bw() +
    facet_wrap(~pheno_cat, nrow = 1, scale = "free_y", space = "free_y") +
    geom_col(aes(fill = value), width = 1, color = "black", linewidth = 0.4) +
    scale_fill_paletteer_d("ggthemes::Tableau_10", name = NULL) +
    labs(y = NULL, x = xlab) +
    scale_x_continuous(expand = expansion(0)) +
    theme(strip.text = element_text(hjust = 0, face = "bold", size = 9,
                                    margin = margin(b = 0.2, unit = "lines")),
          strip.background = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          legend.text = element_text(size = 9, margin = margin(l = 2)),
          legend.key.spacing.y = unit(0.25, "lines"),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.2, "lines"))
}

# --- MAGMA vs. MR gene summary by trait ---

(p_gwas_n <- df_gene_summary_by_pheno %>%
    filter(annotation == "gwas") %>%
    mkplot_gene_summary(x_col = n_genes,
                        xlab = "Number of unique genes") +
    theme(plot.margin = margin(),
          legend.position = "none")
)

(p_gwas_prop <- df_gene_summary_by_pheno %>%
    filter(annotation == "gwas") %>%
    mkplot_gene_summary(x_col = prop_genes,
                        xlab = "Proportion of genes") +
    scale_x_continuous(labels = percent, expand = expansion(0)) +
    scale_fill_paletteer_d("ggthemes::Tableau_10",
      breaks = c("TenK10K MR only", "TenK10K MR & MAGMA", "MAGMA only"),
      name = NULL
    ) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(color = NA),
          axis.ticks.y = element_line(color = "gray"),
          legend.direction = "horizontal",
          legend.position = "inside",
          legend.position.inside = c(1, 1),
          legend.justification = c(1, 0),
          plot.margin = margin(r = 0.5, unit = "lines"),
          axis.ticks.length.y = unit(1, "lines"))
)

p_gwas_mr <- p_gwas_n + p_gwas_prop +
  plot_layout(nrow = 1, guides = "keep") &
  layer_scales(p_gwas_n)$y

# Supplementary Fig. 10
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig-s10.png",
       p_gwas_mr, scaling = 1.05,
       width = 8, height = 12, device = agg_png)
# cairo_pdf has no `scaling` argument; scale = 1 / 1.05 keeps the same text-to-canvas ratio
ggsave("figures/fig-s10.pdf",
       p_gwas_mr, scale = 1 / 1.05,
       width = 8, height = 12, device = cairo_pdf)

# --- eQTLgen MR vs. TenK10K MR gene summary by trait ---

(p_eqtlgen_n <- df_gene_summary_by_pheno %>%
    filter(annotation == "eqtlgen") %>%
    mkplot_gene_summary(x_col = n_genes,
                        xlab = "Number of unique genes") +
    theme(plot.margin = margin(),
          legend.position = "none")
)

(p_eqtlgen_prop <- df_gene_summary_by_pheno %>%
    filter(annotation == "eqtlgen") %>%
    mkplot_gene_summary(x_col = prop_genes,
                        xlab = "Proportion of genes") +
    scale_x_continuous(labels = percent, expand = expansion(0)) +
    scale_fill_paletteer_d("ggthemes::Tableau_10",
      breaks = c("TenK10K MR only", "TenK10K MR & eQTLgen MR", "eQTLgen MR only"),
      name = NULL
    ) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(color = NA),
          axis.ticks.y = element_line(color = "gray"),
          legend.direction = "horizontal",
          legend.position = "inside",
          legend.position.inside = c(1, 1),
          legend.justification = c(1, 0),
          plot.margin = margin(r = 0.5, unit = "lines"),
          axis.ticks.length.y = unit(1, "lines"))
)

(p_eqtlgen <- p_eqtlgen_n + p_eqtlgen_prop +
    plot_layout(nrow = 1, guides = "keep") &
    layer_scales(p_eqtlgen_n)$y
)

# Supplementary Fig. 11
ggsave("figures/fig-s11.png",
       p_eqtlgen, scaling = 1.05,
       width = 8, height = 12, device = agg_png)
ggsave("figures/fig-s11.pdf",
       p_eqtlgen, scale = 1 / 1.05,
       width = 8, height = 12, device = cairo_pdf)
