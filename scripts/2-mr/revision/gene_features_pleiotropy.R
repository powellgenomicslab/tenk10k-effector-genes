library(data.table)
library(tidyverse)
library(arrow)
library(patchwork)
library(ragg)
library(ggrepel)
library(readxl)
library(paletteer)
library(scales)
library(geomtextpath)

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v3.parquet.gz")
# calculate N phenotype / phenotype category

df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")
df_cell_map <- fread("resources/metadata/cell_map.tsv")
INPUT <- list(
  gen_cor = "results/aggregate/tenk10k_phase1.gen_cor.ldak.tsv",
  mr_rrho = "results/rrho/spearman_corr_all_trait_combos_strictmr.tsv"
)
# Create phenotype correlation matrix (same as in calc_indep_clust.R)
# similar analysis using RRHO correlation estimates
df_rrho <- read_tsv(INPUT[["mr_rrho"]]) %>% 
  select(trait1 = p1, trait2 = p2, value = spearman_corr)

# Build square symmetric matrix with 1s on the diagonal
rrho_traits <- union(df_rrho$trait1, df_rrho$trait2)
mat_cor <- bind_rows(
    df_rrho,
    rename(df_rrho, trait1 = trait2, trait2 = trait1),
    tibble(trait1 = rrho_traits, trait2 = rrho_traits, value = 1)
  ) %>%
  pivot_wider(names_from = trait2, values_from = value) %>%
  column_to_rownames("trait1") %>%
  data.matrix() %>%
  .[rrho_traits, rrho_traits]

gene_summary <- df_msmr_tenk10k[mr == TRUE, .(
  n_phenotypes = length(unique(pheno_label)),
  n_disease_phenotypes = length(unique(pheno_label[supercategory == "disease"])),
  n_pheno_cat = length(unique(pheno_cat)),
  n_celltypes = length(unique(biosample)),
  n_major_celltypes = length(unique(major_cell_type)),
  phenotypes = list(unique(pheno_label)),
  celltypes = list(unique(biosample)),
  min_phet_ivw = min(phet_ivw, na.rm = TRUE) |> na_if(Inf),
  min_psigmay_mrlink2 = min(psigmay_mrlink2, na.rm = TRUE) |> na_if(Inf),
  min_p_HEIDI = min(p_HEIDI, na.rm = TRUE) |> na_if(Inf)
), by = probeID]


# Function using eigenvalue decomposition
count_effective_tests <- function(pheno_vec, cor_matrix) {
  if (length(pheno_vec) <= 1) return(length(pheno_vec))
  
  pheno_in_matrix <- intersect(pheno_vec, rownames(cor_matrix))
  if (length(pheno_in_matrix) <= 1) return(length(pheno_in_matrix))
  
  sub_cor <- cor_matrix[pheno_in_matrix, pheno_in_matrix]
  
  eigenvalues <- eigen(sub_cor, symmetric = TRUE, only.values = TRUE)$values
  
  # Li and Ji (2005) method
  f <- function(x) {
    if (x <= 0) return(0)
    indicator <- as.numeric(x >= 1)
    fractional <- x - floor(x)
    return(indicator + fractional)
  }
  
  meff <- sum(sapply(eigenvalues, f))
  return(meff)
}

# Apply to gene summary
gene_summary[, n_effective_phenotypes := sapply(phenotypes, function(p) {
  count_effective_tests(p, mat_cor)
})]

# overall N effective phenotypes
phenos <- df_msmr_tenk10k[mr == TRUE, unique(pheno_label)]
overall_meff <- count_effective_tests(phenos, mat_cor)


# Constraints
# constraint_metrics <- fread("resources/misc/gnomad.v4.1.constraint_metrics.tsv") |> 
#   filter(str_detect(gene_id, "^ENSG") & canonical == TRUE) |> 
#     select(gene, gene_id, lof.pLI, lof.oe, lof.oe_ci.upper, constraint_flags)
# gene_summary[constraint_metrics, loeuf := i.lof.oe_ci.upper, on = c("probeID" = "gene_id")]
  
# annotate
gene_summary[df_gene_annot, gene := i.hgnc_symbol, on = c("probeID" = "ensembl_gene_id")]
# Extended MHC (xMHC): chr6:25,726,063-33,400,644 GRCh38 (HIST1H2AA to RPL12P1)
# Horton et al. (2004) Nat Rev Genet 5:889-899; GRCh38 coords via https://cloufield.github.io/gwaslab/HLA/
mhc_genes <- df_gene_annot[chr == 6 & start < 33400644 & end > 25726063, hgnc_symbol]
# gene_summary[, gene_label := ifelse(frank(-n_effective_phenotypes) < 20, gene, ""), by = high_ld]

gene_summary[,.N, by = n_celltypes]

(p_top_pleiotropic_genes <- 
  gene_summary |> 
  filter(n_celltypes == 1) |> 
  ggplot(aes(x = fct_reorder(gene, n_effective_phenotypes),
         y = n_effective_phenotypes)) +
  geom_point(aes(fill = factor(n_pheno_cat)), shape = "circle filled", alpha = 1) +
  theme_bw(base_family = "Helvetica") +
  geom_text_repel(
    aes(label = ifelse(frank(-n_effective_phenotypes) < 30, gene, "")),
    size = 8/.pt, direction = "y", force = 4,
    xlim = c(NA, 1000), ylim = c(NA, NA),
    max.overlaps = Inf, fontface = "italic", segment.size = 0.1,
  ) +
  scale_fill_viridis_d(name = "Number of\nPhenotype Categories", direction = -1) +
  theme(panel.grid = element_blank()) +
  labs(x = "Gene", y = "Number of Effective Phenotypes") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
      legend.position = "inside",
      panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    legend.position.inside = c(0, 1),
  legend.justification.inside = c(0, 1),
  legend.background = element_blank(),
  strip.background = element_blank()
))


# Top 10 pleiotropic non-MHC genes for labelling (panel B)
label_genes_b <- gene_summary[n_celltypes == 1 & !gene %in% mhc_genes][order(-n_effective_phenotypes)][1:10, gene]
# quantify average pleiotropy (for text)
gene_summary[,median(n_effective_phenotypes)]


(p_cell_phenotypes <- gene_summary |>
  # filter(min_p_HEIDI > 0.05 | min_phet_ivw > 0.05 | min_psigmay_mrlink2 > 0.05) |>
  # mutate(decile = ntile(n_effective_phenotypes, 10)) |>
  ggplot(aes(x = n_celltypes, y = n_effective_phenotypes, group = n_celltypes)) +
  geom_boxplot(fill = "gray", outliers = FALSE, data = ~.x[n_celltypes == 1]) +
  geom_boxplot(fill = "gray", outliers = TRUE, data = ~.x[n_celltypes > 1]) +
  # Unlabelled outlier points (behind)
  geom_point(data = ~.x[n_celltypes == 1 & !gene %in% label_genes_b][
      which(n_effective_phenotypes %in% boxplot(n_effective_phenotypes, plot=FALSE)$out)],
    aes(fill = unlist(celltypes)), shape = "circle filled",
     alpha = 0.8, size = 2) +
  # Labelled points plotted on top for accurate color
  geom_point(data = ~.x[n_celltypes == 1 & gene %in% label_genes_b],
    aes(fill = unlist(celltypes)), shape = "circle filled",
     alpha = 0.8, size = 2) +
  # Label top pleiotropic non-MHC region genes (1 cell type)
  geom_text_repel(
    data = ~.x[n_celltypes == 1 & gene %in% label_genes_b],
    aes(label = gene),
    size = 6/.pt, direction = "y", force = 3,
    xlim = c(NA, -2), ylim = c(NA, NA), hjust = 1,
    nudge_x = -4, segment.color = "grey50", point.padding = 0.3,
    box.padding = unit(0.25, "lines"),
    max.overlaps = Inf, fontface = "italic", segment.size = 0.1,
  ) +
  theme_bw(base_family = "Helvetica") +
  scale_fill_manual(
    name = "Cell Type",
    values = deframe(df_cell_map[,.(wg2_scpred_prediction, color)]),
    labels = deframe(df_cell_map[,.(wg2_scpred_prediction, cell_type)])) +
  guides(fill = guide_legend(ncol = 3, override.aes = list(size = 3))) +
  scale_x_continuous(breaks = c(0, 10, 20), expand = expansion(add = c(13, 0.5))) +
  scale_y_continuous(limits = c(0, NA), breaks = function(x) c(0, labeling::extended(x[1], x[2], 5)), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of Cell Types", y = "Number of\nIndependent Phenotypes") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        legend.key.spacing.x = unit(0.3, "lines"),
        legend.key.spacing.y = unit(0.2, "lines"),
        legend.margin = margin(0, 0, 0, 0)))

ggsave("figures/mr_overview/top_pleiotropic_genes.png",
      p_top_pleiotropic_genes, width = 9, height = 5, device = ragg::agg_png)
ggsave("figures/mr_overview/n_phenotypes_by_n_celltypes_filtered_pleio.png",
      p_cell_phenotypes, width = 9, height = 5, device = ragg::agg_png)

saveRDS(p_top_pleiotropic_genes, "figures/mr_overview/top_pleiotropic_genes.rds")
saveRDS(p_cell_phenotypes, "figures/mr_overview/n_phenotypes_by_n_celltypes.rds")


# ------------------------------------------------------------------------------------------------
# Plot n_celltypes vs heterogeneity 
# ------------------------------------------------------------------------------------------------

# gene_summary[, gene_heterogeneity := if_else(min_p_HEIDI > 0.05 | min_phet_ivw > 0.05 | min_psigmay_mrlink2 > 0.05, "Not heterogeneous", "Heterogeneous")]
# gene_summary[, gene_heterogeneity := if_else(is.na(gene_heterogeneity), "NA", gene_heterogeneity)]

# # Compute tally with proportion of heterogeneous genes per n_celltypes
# df_het_tally <- gene_summary[, .N, by = .(n_celltypes, gene_heterogeneity)] |>
#   as_tibble() |>
#   complete(n_celltypes, gene_heterogeneity, fill = list(N = 0)) |>
#   mutate(sum_n = sum(N), prop = N / sum_n, .by = n_celltypes) |>
#   mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))

# (p_n_celltypes_vs_heterogeneity <-
#   ggplot(df_het_tally, aes(x = n_celltypes, y = N)) +
#   geom_col(aes(fill = gene_heterogeneity), width = 1, linewidth = 0.5,
#            position = "stack", color = "black") +
#   scale_fill_manual(values = c("Heterogeneous" = "#4C6C94FF", "Not heterogeneous" = "#A4ABB0FF", "NA" = "grey80"),
#                     name = NULL) +
#   geom_line(aes(y = prop_scaled),
#             data = ~filter(.x, gene_heterogeneity == "Heterogeneous")) +
#   geom_textline(aes(label = "% Heterogeneous", y = prop_scaled),
#                 data = ~filter(.x, gene_heterogeneity == "Heterogeneous"),
#                 size = 9/.pt, color = "red3",
#                 vjust = -1, hjust = 0.4) +
#   geom_point(aes(y = prop_scaled),
#              data = ~filter(.x, gene_heterogeneity == "Heterogeneous"),
#              color = "red3") +
#   theme_classic(base_family = "Helvetica") +
#   labs(x = "Number of Cell Types", y = "Number of MR genes") +
#   scale_x_continuous(breaks = seq(1, 28, 4),
#                      expand = expansion(add = 0.25)) +
#   scale_y_continuous(labels = label_comma(),
#                      n.breaks = 10,
#                      expand = expansion(mult = c(0, 0.05)),
#                      sec.axis = sec_axis(
#                        ~rescale(., to = c(0, 1)),
#                        name = NULL,
#                        labels = label_percent(accuracy = 1)),
#                      guide = guide_axis(cap = "none")) +
#   coord_cartesian(clip = "off") +
#   theme(legend.position.inside = c(0.95, 1),
#         legend.position = "inside",
#         axis.title = element_text(size = 9),
#         legend.text = element_text(size = 9, margin = margin(l = 2)),
#         legend.justification = c(1, 1),
#         legend.key.size = unit(0.75, "lines"),
#         legend.key.spacing.y = unit(0.5, "lines"),
#         axis.text = element_text(size = 8),
#         axis.line.y.right = element_line(color = "red3"),
#         axis.text.y.right = element_text(color = "red3"),
#         axis.ticks.y.right = element_line(color = "red3"),
#         axis.title.y.right = element_text(color = "red3"),
#         panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
#   ))

# p_n_celltypes_vs_heterogeneity %>%
#   ggsave(filename = "figures/mr_overview/n_celltypes_vs_heterogeneity.png", width = 9, height = 5, device = ragg::agg_png)

# # Compute tally with proportion of heterogeneous genes per n_effective_phenotypes
# df_het_tally_pheno <- gene_summary[, .(n_effective_phenotypes_int = round(n_effective_phenotypes))
#                                    ][, cbind(gene_summary, .SD)] |>
#   _[, .N, by = .(n_effective_phenotypes_int, gene_heterogeneity)] |>
#   as_tibble() |>
#   complete(n_effective_phenotypes_int, gene_heterogeneity, fill = list(N = 0)) |>
#   mutate(sum_n = sum(N), prop = N / sum_n, .by = n_effective_phenotypes_int) |>
#   mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))

# (p_n_phenotypes_vs_heterogeneity <-
#   ggplot(df_het_tally_pheno, aes(x = n_effective_phenotypes_int, y = N)) +
#   geom_col(aes(fill = gene_heterogeneity), width = 1, linewidth = 0.5,
#            position = "stack", color = "black") +
#   scale_fill_manual(values = c("Heterogeneous" = "#4C6C94FF", "Not heterogeneous" = "#A4ABB0FF", "NA" = "grey80"),
#                     name = NULL) +
#   geom_line(aes(y = prop_scaled),
#             data = ~filter(.x, gene_heterogeneity == "Heterogeneous")) +
#   geom_textline(aes(label = "% Heterogeneous", y = prop_scaled),
#                 data = ~filter(.x, gene_heterogeneity == "Heterogeneous"),
#                 size = 9/.pt, color = "red3",
#                 vjust = -1, hjust = 0.4) +
#   geom_point(aes(y = prop_scaled),
#              data = ~filter(.x, gene_heterogeneity == "Heterogeneous"),
#              color = "red3") +
#   theme_classic(base_family = "Helvetica") +
#   labs(x = "Number of Effective Phenotypes", y = "Number of MR genes") +
#   scale_x_continuous(expand = expansion(add = 0.25)) +
#   scale_y_continuous(labels = label_comma(),
#                      n.breaks = 10,
#                      expand = expansion(mult = c(0, 0.05)),
#                      sec.axis = sec_axis(
#                        ~rescale(., to = c(0, 1)),
#                        name = NULL,
#                        labels = label_percent(accuracy = 1)),
#                      guide = guide_axis(cap = "none")) +
#   coord_cartesian(clip = "off") +
#   theme(legend.position.inside = c(0.95, 1),
#         legend.position = "inside",
#         axis.title = element_text(size = 9),
#         legend.text = element_text(size = 9, margin = margin(l = 2)),
#         legend.justification = c(1, 1),
#         legend.key.size = unit(0.75, "lines"),
#         legend.key.spacing.y = unit(0.5, "lines"),
#         axis.text = element_text(size = 8),
#         axis.line.y.right = element_line(color = "red3"),
#         axis.text.y.right = element_text(color = "red3"),
#         axis.ticks.y.right = element_line(color = "red3"),
#         axis.title.y.right = element_text(color = "red3"),
#         panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
#   ))

# p_n_phenotypes_vs_heterogeneity %>%
#   ggsave(filename = "figures/mr_overview/n_phenotypes_vs_heterogeneity.png", width = 9, height = 5, device = ragg::agg_png)

# (p_heterogeneity_combined <- p_n_celltypes_vs_heterogeneity + p_n_phenotypes_vs_heterogeneity +
#   plot_layout(ncol = 2, guides = "collect") +
#   plot_annotation(tag_levels = "a") &
#   theme(plot.tag = element_text(size = 14, face = "bold"),
#         legend.position = "bottom"))

# ggsave("figures/mr_overview/heterogeneity_combined.png",
#        p_heterogeneity_combined, width = 14, height = 5, device = ragg::agg_png)

# ------------------------------------------------------------------------------------------------
# Plot the cell type specific pleiotropic genes
# ------------------------------------------------------------------------------------------------

# get the pleitropic genes that are also highly cell type specific 
top_pleiotropic_genes <- gene_summary |> 
  filter(min_p_HEIDI > 0.05 | min_phet_ivw > 0.05 | min_psigmay_mrlink2 > 0.05) |>
  filter(n_major_celltypes == 1 & n_celltypes > 1) %>% 
  slice_max(n_disease_phenotypes, n=8)

selected_genes <- c("ERBB2", "MYRF")

plot_data <- df_msmr_tenk10k  %>% 
    # filter(Gene %in% unique(top_pleiotropic_genes$gene))  %>% 
    # filter(Gene %in% selected_genes, supercategory == "disease") %>%
    filter(Gene %in% selected_genes) %>%
    mutate(signed_neg_log10_p = -log10(p_SMR_multi) * sign(b_SMR)) %>% 
    group_by(Gene) %>%
    mutate(
      x = frank(-signed_neg_log10_p, ties.method = "average")
    ) %>% 
    ungroup() %>% 
    mutate(
      # Significant = case_when(mr_plus & mr ~ "MR+", mr ~ "MR", !mr ~ "Not significant"),
      Significant = if_else(is.na(max_evidence), "Not significant", max_evidence),
      `Phenotype category` = pheno_cat,
    ) # note only works if gene is only tested in one cell type

label_data <- plot_data %>%
  filter(mr & supercategory == "disease") %>%
  mutate(mid_x = max(x) / 2, .by = Gene)

# color palette 
# color_df <- fread("resources/misc/trait_cat_col.tsv")
# disease_category_cols <- setNames(color_df$value[which(color_df$name %in% unique(plot_data$pheno_cat))], color_df$name[which(color_df$name %in% unique(plot_data$pheno_cat))])
# # bio_category_cols <- setNames(color_df$value[which(color_df$name %in% unique(bio_trait_meta$pheno_cat))], color_df$name[which(color_df$name %in% unique(bio_trait_meta$pheno_cat))])

cat_order <- read_xlsx("resources/metadata/trait_metadata_curated.xlsx",
  sheet = "trait_category_order"
) %>%
  pull(cat_order)

trait_cat_col <- paletteer_d("ggthemes::Tableau_10", 10) %>%
  set_names(cat_order)
disease_category_cols <- trait_cat_col[names(trait_cat_col) %in% unique(plot_data$`Phenotype category`)]

# Todo:
# [x] only label points for biological phenotypes 
# [] change the shape scale to be more intuitive. 

# Label data for inside-panel strip text (top-right corner)
facet_label_data <- plot_data %>%
  distinct(Gene, cell_type) %>%
  mutate(label = paste0("\"", cell_type, ":\" ~ italic(\"", Gene, "\")"))

p_scatter_pleiotropic <- plot_data %>%
    ggplot(aes(shape = Significant, y = signed_neg_log10_p, x = x, color=`Phenotype category`, fill = `Phenotype category`)) +
    facet_wrap(vars(Gene, cell_type), scales = "free", nrow = 1) +
    geom_point(data = plot_data %>% filter(!mr), color = "lightgrey", fill = "lightgrey", size = rel(2),
      show.legend = c(shape = TRUE, color = FALSE, fill = FALSE)) +
    geom_point(data = plot_data %>% filter(mr), size = rel(2)) +
    geom_text_repel(data = label_data,
      aes(label = pheno_label),
      nudge_x = label_data$mid_x - label_data$x,
      size = rel(3),
      direction = "y",
      min.segment.length = 0,
      hjust = 0.5,
      segment.color = "grey50",
      point.padding = 0.5,
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.size = 0.1,
      force = 4,
      show.legend = FALSE
    ) +
    geom_text(data = facet_label_data,
      aes(x = Inf, y = Inf, label = label),
      parse = TRUE, inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.3,
      size = 8 * 1.3 / .pt, family = "Helvetica") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    guides(
      y = guide_axis(cap = "none"),
      colour = guide_legend(ncol = 1),
      fill = guide_legend(ncol = 1),
      shape = guide_legend(ncol = 1, override.aes = list(fill = "black"))
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 8, base_family = "Helvetica") +
    theme(
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    # centre on zero
    labs(y = expression(sign(beta)%*%-log[10](p[SMR_multi]))) +
    scale_color_manual(values = disease_category_cols) +
    scale_fill_manual(values = disease_category_cols) +
    # MR : circle, MR+ : diamond, not significant : cross
    scale_shape_manual(
      values = c(
        "Not significant" = 3,
        "mr" = 1,
        "mr_coloc" = 2,
        "mr_sens" =  21,
        "mr_sens_coloc" = 24)) +
    theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "vertical",
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        strip.clip = "off")

p_scatter_pleiotropic %>% 
  ggsave(filename = "fig_pub/fig_1_pleiotropic_genes_scatter.png", width = 2 * length(unique(plot_data$Gene)), height = 8, limitsize = FALSE)

p_scatter_pleiotropic %>%
  ggsave(filename = "fig_pub/fig_1_pleiotropic_genes_scatter.pdf", width = 2 * length(unique(plot_data$Gene)), height = 8, limitsize = FALSE)

# Apply combined-figure theme adjustments before saving RDS
p_scatter_pleiotropic <- p_scatter_pleiotropic +
  guides(
    colour = guide_legend(ncol = 1),
    shape  = guide_legend(ncol = 1),
    size   = guide_legend(ncol = 1)
  ) +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(0.75, "lines"),
    legend.key.spacing.y = unit(0.25, "lines"),
    plot.title = element_blank(),
  )

p_scatter_pleiotropic %>% saveRDS("fig_pub/fig_1_pleiotropic_genes_scatter.rds")

# ------------------------------------------------------------------------------------------------
# Volcano plot: beta_SMR vs -log10(p_SMR_multi)
# ------------------------------------------------------------------------------------------------

# p_volcano <- plot_data %>%
#     ggplot(aes(x = b_SMR, y = -log10(p_SMR_multi), color = `Phenotype category`, shape = Significant)) +
#     geom_point(data = plot_data %>% filter(!mr), color = "lightgrey", size = rel(2)) +
#     geom_point(data = plot_data %>% filter(mr), size = rel(2)) +
#     geom_text_repel(
#       data = plot_data %>% filter(mr),
#       aes(label = pheno_label), size = rel(2),
#       min.segment.length = 0,
#       segment.color = "grey50",
#       point.padding = 0.3,
#       max.overlaps = 20,
#       segment.size = 0.1,
#       show.legend = FALSE
#     ) +
#     facet_wrap(vars(Gene, cell_type), scales = "free", strip.position = "bottom", nrow = 1,
#       labeller = label_bquote(.(as.character(cell_type)) * ":" ~ italic(.(Gene)))) +
#     scale_color_manual(values = disease_category_cols) +
#     scale_x_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
#     theme_classic() +
#     labs(x = expression(beta[SMR]), y = expression(-log[10](p[SMR_multi]))) +
#     theme(
#       strip.placement = "outside",
#       strip.background = element_blank(),
#       strip.clip = "off",
#       panel.grid.major = element_line(color = "gray90")
#     )

# p_volcano %>%
#   ggsave(filename = "fig_pub/volcano_pleiotropic_genes.png",
#     width = 6 * length(unique(plot_data$Gene)), height = 6, limitsize = FALSE)

# p_volcano %>%
#   ggsave(filename = "fig_pub/volcano_pleiotropic_genes.pdf",
#     width = 6 * length(unique(plot_data$Gene)), height = 6, limitsize = FALSE)

# p_volcano %>% saveRDS("fig_pub/volcano_pleiotropic_genes.rds")

# # ------------------------------------------------------------------------------------------------
# # Lollipop plot: b_SMR with point size = -log10(p_SMR_multi)
# # ------------------------------------------------------------------------------------------------

# lollipop_data <- plot_data %>%
#   mutate(
#     neg_log10_p = -log10(p_SMR_multi),
#     signed_neg_log10_p = neg_log10_p * sign(b_SMR),
#     abs_beta = abs(b_SMR)
#   ) %>%
#   arrange(Gene, cell_type, signed_neg_log10_p) %>%
#   group_by(Gene, cell_type) %>%
#   mutate(rank = row_number()) %>%
#   ungroup()

# p_lollipop <- lollipop_data %>%
#     ggplot(aes(x = rank, y = signed_neg_log10_p)) +
#     geom_segment(
#       aes(xend = rank, y = 0, yend = signed_neg_log10_p),
#       color = "grey70", linewidth = 0.4
#     ) +
#     geom_point(
#       aes(size = abs_beta),
#       data = lollipop_data %>% filter(!mr),
#       color = "lightgrey"
#     ) +
#     geom_point(
#       aes(color = `Phenotype category`, size = abs_beta),
#       data = lollipop_data %>% filter(mr)
#     ) +
#     geom_text_repel(
#       data = lollipop_data %>% filter(mr),
#       aes(label = pheno_label, color = `Phenotype category`), size = rel(2),
#       direction = "y",
#       nudge_x = 10,
#       min.segment.length = 0,
#       segment.color = "grey50",
#       point.padding = 0.3,
#       max.overlaps = Inf,
#       segment.size = 0.1,
#       show.legend = FALSE
#     ) +
#     scale_color_manual(values = disease_category_cols) +
#     scale_size_continuous(name = expression("|" * beta[SMR] * "|"), range = c(0.5, 4)) +
#     coord_cartesian(clip = "off") +
#     facet_wrap(vars(Gene, cell_type), scales = "free", strip.position = "bottom", ncol = 1,
#       labeller = label_bquote(.(as.character(cell_type)) * ":" ~ italic(.(Gene)))) +
#     theme_classic() +
#     labs(x = NULL, y = expression(sign(beta[SMR]) %*% -log[10](p[SMR_multi]))) +
#     theme(
#       strip.placement = "outside",
#       strip.background = element_blank(),
#       strip.clip = "off",
#       panel.grid.major.y = element_line(color = "gray90"),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.line.x = element_blank()
#     )

# p_lollipop %>%
#   ggsave(filename = "fig_pub/lollipop_pleiotropic_genes.png",
#     width = 4 * length(unique(plot_data$Gene)), height = 3 * length(unique(plot_data$Gene)), limitsize = FALSE)

# p_lollipop %>%
#   ggsave(filename = "fig_pub/lollipop_pleiotropic_genes.pdf",
#     width = 4 * length(unique(plot_data$Gene)), height = 3 * length(unique(plot_data$Gene)), limitsize = FALSE)

# p_lollipop %>% saveRDS("fig_pub/lollipop_pleiotropic_genes.rds")