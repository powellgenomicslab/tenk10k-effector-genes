# figures: MR vs. MAGMA correlation by trait, and MR test count vs. hit count by cell type
#
# Recovered from scripts/mr_overview/overview_supplement.R (deleted in commit e137f97,
# 2026-05-12). Two figures:
#   - fig-s9 (smr_gwas_correlation): number of MAGMA genes vs. number of MR genes, per trait
#   - fig-s3 (smr_celltype_scatter): number of MR tests vs. number of significant MR hits, per cell type
#
# The df_summary_gene_pheno build below matches scripts/table/tab_mr_overlap.R's
# mr_gwas_overlap table build, kept in sync so the figure and the table agree.

source("scripts/0-preprocess/preprocess_results.R")

library(scales)
library(ragg)
library(paletteer)

df_msmr <- df_msmr_tenk10k[mr == TRUE]

# --- MAGMA vs. MR gene count correlation, per trait ---

df_summary_gene_pheno <- df_msmr_tenk10k %>%
  group_by(pheno_label, pheno_cat, probeID) %>%
  summarise(mr = max(mr), gwas = max(magma_gene), .groups = "drop_last") %>%
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
    jaccard_mr_gwas = sum(mr & gwas) / sum(mr | gwas),
    .groups = "drop"
  ) %>%
  setDT()

rho_gwas_mr <- cor.test(df_summary_gene_pheno$n_gwas, df_summary_gene_pheno$n_mr)
print(rho_gwas_mr)

(p_cor_gwas_mr <- ggplot(df_summary_gene_pheno,
       aes(x = n_gwas, y = n_mr)) +
  theme_bw() +
  geom_point(aes(color = pheno_cat), size = 3, alpha = 0.8) +
  geom_smooth(color = "gray20", linewidth = 0.5) +
  paletteer::scale_color_paletteer_d("ggthemes::Tableau_10",
                                      name = "Phenotype category") +
  labs(x = "Number of MAGMA genes",
       y = "Number of MR genes") +
  scale_x_continuous(expand = 0.01) +
  scale_y_continuous(expand = 0.01) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linewidth = 0.5, linetype = "dashed") +
  coord_equal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0, 1)
        )
)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave(p_cor_gwas_mr, filename = "figures/fig-s9.png",
       width = 6, height = 8, device = agg_png, scaling = 1.2)
ggsave(p_cor_gwas_mr, filename = "figures/fig-s9.pdf",
       width = 6 / 1.2, height = 8 / 1.2, device = cairo_pdf)

# --- MR test count vs. significant MR hit count, per cell type ---

df_cell_summary <- df_msmr_tenk10k[,
  .(n_genes = length(unique(probeID)),
    n_mr_genes = .SD[mr == TRUE, length(unique(probeID))],
    n_sig = sum(mr),
    n_tested = .N,
    prop_sig = sum(mr) / .N
    ),
  by = list(cell_type, major_cell_type)]

(p_scatter_cell <- df_cell_summary %>%
  ggplot(aes(x = n_sig, y = n_tested)) +
  theme_bw() +
  geom_smooth(color = "gray20", linewidth = 0.5) +
  geom_point(aes(color = cell_type)) +
  scale_color_manual(
    values = deframe(df_cell_map[,.(cell_type, color)]),
    name = NULL) +
  labs(x = "Number of significant MR associations",
       y = "Number of MR tests") +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 9))
)

ggsave("figures/fig-s3.png",
       p_scatter_cell, scaling = 1.5,
       width = 10, height = 6, device = agg_png)
ggsave("figures/fig-s3.pdf",
       p_scatter_cell,
       width = 10 / 1.5, height = 6 / 1.5, device = cairo_pdf)
