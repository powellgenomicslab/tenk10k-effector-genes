library(data.table)
library(tidyverse)
library(arrow)
library(patchwork)
library(ragg)
library(ggrepel)
library(paletteer)

dir.create("figures/mr_overview", recursive = TRUE, showWarnings = FALSE)

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v5.parquet.gz")

df_gene_annot <- fread("metadata/gencode.v44.gene_type.tsv.gz")
df_cell_map <- fread("metadata/cell.tsv")
INPUT <- list(
  mr_spearman = "results/rrho/spearman_corr_all_trait_combos_strictmr.tsv"
)
# Create phenotype correlation matrix
# using pairwise spearman correlations of MR results between traits
df_spearman <- read_tsv(INPUT[["mr_spearman"]]) %>%
  select(trait1 = p1, trait2 = p2, value = spearman_corr)

# Build square symmetric matrix with 1s on the diagonal
spearman_traits <- union(df_spearman$trait1, df_spearman$trait2)
mat_cor <- bind_rows(
    df_spearman,
    rename(df_spearman, trait1 = trait2, trait2 = trait1),
    tibble(trait1 = spearman_traits, trait2 = spearman_traits, value = 1)
  ) %>%
  pivot_wider(names_from = trait2, values_from = value) %>%
  column_to_rownames("trait1") %>%
  data.matrix() %>%
  .[spearman_traits, spearman_traits]

gene_summary <- df_msmr_tenk10k[mr == TRUE, .(
  n_phenotypes = length(unique(pheno_label)),
  n_celltypes = length(unique(biosample)),
  phenotypes = list(unique(pheno_label)),
  celltypes = list(unique(biosample))
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

# annotate
gene_summary[df_gene_annot, gene := i.hgnc_symbol, on = c("probeID" = "ensembl_gene_id")]
# Extended MHC (xMHC): chr6:25,726,063-33,400,644 GRCh38 (HIST1H2AA to RPL12P1)
# Horton et al. (2004) Nat Rev Genet 5:889-899; GRCh38 coords via https://cloufield.github.io/gwaslab/HLA/
mhc_genes <- df_gene_annot[chr == 6 & start < 33400644 & end > 25726063, hgnc_symbol]

gene_summary[,.N, by = n_celltypes]

# Top 10 pleiotropic non-MHC genes for labelling (panel B)
label_genes_b <- gene_summary[n_celltypes == 1 & !gene %in% mhc_genes][order(-n_effective_phenotypes)][1:10, gene]
# top of the n_celltypes == 1 genes: value, cell type and label / MHC status
print(gene_summary[n_celltypes == 1][order(-n_effective_phenotypes)][1:20, .(
  gene,
  n_effective_phenotypes,
  cell_type = unlist(celltypes),
  labelled = gene %in% label_genes_b,
  mhc = gene %in% mhc_genes
)])
# quantify average pleiotropy (for text)
gene_summary[,median(n_effective_phenotypes)]

# ------------------------------------------------------------------------------------------------
# Panel b: genes in 1 cell type by cell type (top, labelled) over genes in 2+ cell types by number of cell types (bottom)
# ------------------------------------------------------------------------------------------------

hla_colours <- c(`TRUE` = "#D55E00", `FALSE` = "#0072B2") # Okabe-Ito vermillion / blue
hla_labels <- c(`TRUE` = "HLA region gene", `FALSE` = "non-HLA region gene")

one_celltype <- gene_summary |>
  filter(n_celltypes == 1) |>
  mutate(cell_type = fct_reorder(unlist(celltypes), n_effective_phenotypes, .desc = TRUE),
         hla = gene %in% mhc_genes)

# boxplot outliers drawn as a separate layer so they can carry the HLA colour
one_celltype_outliers <- one_celltype |>
  filter(n_effective_phenotypes %in% boxplot.stats(n_effective_phenotypes)$out, .by = cell_type)

# label positions hard-coded (row 1 = lower, row 2 = upper; x_label = left edge of the text in x-axis positions);
# positions are precomputed columns so the saved plot has no free variables. Genes missing from this table are not labelled.
label_pos_one <- tribble(
  ~gene,             ~row, ~x_label, ~x_segment, # x_segment: where the leader leaves the bottom of the text
  "INKA1",              1,      1,       2.75,
  "RPL6",               1,    4.8,        6.2,
  "ERBB2",              1,    8.3,         10,
  "PGAP3",              1,   15.8,       17.5,
  "P4HA2",              1,     20,       21.75,
  "HECTD4",             2,      1,        4.6,
  "MAPT",               2,    5.7,        8.1,
  "MYRF",               2,   10.5,       12.9,
  "GPBAR1",             2,     14,       16.1,
  "ENSG00000254936",    2,   18.3,       26.8
)
y_top_one <- max(one_celltype$n_effective_phenotypes)
label_data_one <- one_celltype |>
  filter(gene %in% label_genes_b) |>
  as_tibble() |>
  inner_join(label_pos_one, by = "gene") |>
  mutate(y_label = y_top_one * c(1.14, 1.38)[row], # row heights relative to the data maximum
         y_segment = y_label - 0.05 * y_top_one)

p_cell_phenotypes_one <- one_celltype |>
  ggplot(aes(x = cell_type, y = n_effective_phenotypes)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(data = one_celltype_outliers, aes(colour = hla), size = 0.6) +
  geom_segment(data = label_data_one,
    aes(x = x_segment, xend = as.integer(cell_type), y = y_segment, yend = n_effective_phenotypes),
    colour = "grey50", linewidth = 0.15) +
  geom_point(data = label_data_one, size = 0.6, colour = hla_colours[["FALSE"]]) +
  geom_text(data = label_data_one, aes(x = x_label, y = y_label, label = gene),
    size = 5/.pt, family = "Helvetica", fontface = "italic", hjust = 0) +
  scale_colour_manual(values = hla_colours, labels = hla_labels, name = NULL) +
  scale_x_discrete(labels = deframe(df_cell_map[,.(wg2_scpred_prediction, cell_type)])) +
  scale_y_continuous(limits = c(0, NA), breaks = function(x) c(0, labeling::extended(x[1], x[2], 5)), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Genes with MR evidence in 1 cell type", x = "Cell type", y = "Number of\nIndependent Phenotypes") +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(axis.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 6, face = "bold"))

multi_celltype <- gene_summary |>
  filter(n_celltypes > 1) |>
  mutate(hla = gene %in% mhc_genes)

multi_celltype_outliers <- multi_celltype |>
  filter(n_effective_phenotypes %in% boxplot.stats(n_effective_phenotypes)$out, .by = n_celltypes)

p_cell_phenotypes_multi <- multi_celltype |>
  ggplot(aes(x = n_celltypes, y = n_effective_phenotypes, group = n_celltypes)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(data = multi_celltype_outliers, aes(colour = hla), size = 0.6) +
  scale_colour_manual(values = hla_colours, labels = hla_labels, name = NULL) +
  scale_x_continuous(breaks = c(2, 10, 20)) +
  scale_y_continuous(limits = c(0, NA), breaks = function(x) c(0, labeling::extended(x[1], x[2], 5)), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Genes with MR evidence in 2+ cell types", x = "Number of cell types", y = "Number of\nIndependent Phenotypes") +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(axis.text = element_text(size = 5),
        plot.title = element_text(size = 6, face = "bold"))

# stacked vertically to fit the Fig 2 panel b slot (2.35 x 3.9 in); identical y axis titles and the HLA legend are merged
p_cell_phenotypes_split <- p_cell_phenotypes_one / p_cell_phenotypes_multi +
  plot_layout(axis_titles = "collect", guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.5, "lines"),
        legend.margin = margin(0, 0, 0, 0))

ggsave("figures/mr_overview/n_phenotypes_by_n_celltypes_split.png",
      p_cell_phenotypes_split, width = 2.35, height = 3.9, dpi = 400, device = ragg::agg_png)

saveRDS(p_cell_phenotypes_split, "figures/mr_overview/n_phenotypes_by_n_celltypes_split.rds")

# ------------------------------------------------------------------------------------------------
# Extended Data Fig 1: the cell type specific pleiotropic genes
# ------------------------------------------------------------------------------------------------

selected_genes <- c("ERBB2", "MYRF")

plot_data <- df_msmr_tenk10k  %>%
    filter(Gene %in% selected_genes) %>%
    mutate(signed_neg_log10_p = -log10(p_SMR_multi) * sign(b_SMR)) %>%
    group_by(Gene) %>%
    mutate(
      x = frank(-signed_neg_log10_p, ties.method = "average")
    ) %>%
    ungroup() %>%
    mutate(
      Significant = if_else(is.na(max_evidence), "Not significant", max_evidence),
      `Phenotype category` = pheno_cat,
    ) # note only works if gene is only tested in one cell type

label_data <- plot_data %>%
  filter(mr & supercategory == "disease") %>%
  mutate(mid_x = max(x) / 2, .by = Gene)

cat_order <- fread("metadata/trait_category.tsv") %>%
  pull(cat_order)

trait_cat_col <- paletteer_d("ggthemes::Tableau_10", 10) %>%
  set_names(cat_order)
disease_category_cols <- trait_cat_col[names(trait_cat_col) %in% unique(plot_data$`Phenotype category`)]

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
      family = "Helvetica",
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
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    labs(y = expression(sign(beta)%*%-log[10](P[MR]))) +
    scale_color_manual(values = disease_category_cols) +
    scale_fill_manual(values = disease_category_cols) +
    scale_shape_manual(
      name = "Evidence",
      labels = c(
        "Not significant" = "No evidence",
        "mr" = "MR only",
        "mr_sens" = "MR + Sensitivity",
        "mr_coloc" = "MR + Coloc",
        "mr_sens_coloc" = "MR + Sensitivity + Coloc"),
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
  ggsave(filename = "figures/mr_overview/fig_1_pleiotropic_genes_scatter.png", width = 2 * length(unique(plot_data$Gene)), height = 8, limitsize = FALSE)

p_scatter_pleiotropic %>%
  ggsave(filename = "figures/mr_overview/fig_1_pleiotropic_genes_scatter.pdf", device = cairo_pdf, width = 2 * length(unique(plot_data$Gene)), height = 8, limitsize = FALSE)

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

p_scatter_pleiotropic %>% saveRDS("figures/mr_overview/fig_1_pleiotropic_genes_scatter.rds")
