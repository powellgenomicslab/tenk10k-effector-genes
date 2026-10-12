library(data.table)
library(tidyverse)
library(arrow)
library(qvalue)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(ragg)

dir.create("figures/v3/supp", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------------------------------------------------------------
# Load data directly from the v5 parquet
# ----------------------------------------------------------------------------------------------------------------------------------

major_cell_type_order <- c(
  "CD4 T", "CD8 T", "Unconventional T", "NK", "Plasma B", "B", "Monocyte", "Dendritic", "HSPC"
)

df_cell_map <- fread("metadata/cell.tsv") |>
  select(-major_cell_type) |>
  rename(major_cell_type = revision_major_cell_type) |>
  mutate(major_cell_type = factor(major_cell_type, levels = major_cell_type_order))

df_trait_map_all <- fread("metadata/trait.tsv")
df_trait_map <- filter(df_trait_map_all, include)
setDT(df_trait_map)
phenotypes <- df_trait_map$trait_id

cat_order <- fread("metadata/trait_category.tsv") %>%
  pull(cat_order)

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v5.parquet.gz")

gene_universe <- unique(df_msmr_tenk10k$probeID)

df_magma_all <- read_parquet("results/aggregate/tenk10k_phase1.magma.gz.parquet")
df_magma <- df_magma_all %>%
  filter(GENE %in% gene_universe, phenotype %in% phenotypes) %>%
  mutate(qval = qvalue(P)$qvalues,
         lfdr = qvalue(P)$lfdr,
         pbh = p.adjust(P, "BH")) %>%
  filter(lfdr < 0.05) %>%
  setDT()

df_msmr <- df_msmr_tenk10k %>% filter(mr)

# calculate the number of genes tested per cell type / trait

# Matrix for the main heatmap (only significant hits (qval_pheno < 0.05, p_heidi > 0.05))
msmr_mat <- df_msmr %>%
  filter(supercategory == "biological") %>%
  count(cell_type, pheno_label, name = "Number of SMR+HEIDI significant genes") %>%
  pivot_wider(names_from = pheno_label, values_from = `Number of SMR+HEIDI significant genes`) %>%
  column_to_rownames("cell_type") %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  as.matrix()
# order rows
msmr_mat <- msmr_mat[df_cell_map$cell_type, ]

# Matrix for side annotation (eGene / SMRsig counts)
egenes_mat <- df_msmr_tenk10k %>%
  filter(supercategory == "biological") %>%
  mutate(sig = fifelse(lfdr_msmr_pheno < 0.05 & p_HEIDI >= 0.05, "MR eGene", "eGene")) %>%
  select(sig, cell_type, phenotype, probeID) %>%
  # Dedupe the genes. If gene is SMR hit in any cell type x trait combo we will consider/plot that as SMRsig
  mutate(sig = factor(sig, levels = c("MR eGene", "eGene"))) %>%
  arrange(sig) %>%
  group_by(
    cell_type
  ) %>%
  filter(!duplicated(probeID)) %>%
  ungroup() %>%
  count(sig, cell_type) %>%
  inner_join( # make it the same order as the main matrix
    tibble(
      cell_type = factor(rownames(msmr_mat), levels = df_cell_map$cell_type),
      order = 1:length(rownames(msmr_mat))
    )
  ) %>%
  arrange(order) %>%
  select(-order) %>%
  pivot_wider(names_from = sig, values_from = n) %>%
  column_to_rownames("cell_type") %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  as.matrix()

# Matrix for side annotation (cell counts)

df_cell_counts <- df_cell_map

# Matrix for the top annotations (number of GWAS genes)
genes_tested_per_trait_magma <- df_magma %>% # only includes the significant MAGMA genes
  count(phenotype, name = "Number of MAGMA significant genes") %>%
  inner_join(df_trait_map %>% select(phenotype = trait_id, pheno_label = label, pheno_cat = cat_rev)) %>%
  inner_join( # make it the same order as the main matrix
    tibble(
      pheno_label = colnames(msmr_mat),
      order = 1:length(colnames(msmr_mat))
    )
  ) %>%
  arrange(order) %>%
  select(-order) %>%
  mutate(pheno_cat = factor(pheno_cat, levels = cat_order))

# Matrix for top annotation: SMR Hit count by trait and cell type specificity
smr_celltype_specificity <- df_msmr %>%
  select(cell_type, pheno_label, probeID) %>%
  group_by(pheno_label, probeID) %>%
  mutate(n_celltypes_smr_gene = length(unique(cell_type))) %>%
  ungroup() %>%
  mutate(`Number of cell types with shared MR gene` = cut(n_celltypes_smr_gene,
                                                          breaks = c(0, 1, 5, 10, 15, 20, 28),
                                                          labels = c("1", "2-5", "6-10", "11-15", "16-20", "21-28"),
                                                          include.lowest = TRUE
  )) %>%
  select(-cell_type) %>%
  distinct() %>%
  count(pheno_label, `Number of cell types with shared MR gene`, name = "Number of MR genes") %>%
  mutate(`Total MR genes in trait` = sum(`Number of MR genes`), .by = "pheno_label") %>%
  mutate(`Proportion of total MR genes` = `Number of MR genes` / `Total MR genes in trait`)

smr_celltype_specificity_mat <- smr_celltype_specificity %>%
  select(pheno_label, `Number of cell types with shared MR gene`, `Proportion of total MR genes`) %>%
  pivot_wider(names_from = `Number of cell types with shared MR gene`, values_from = `Proportion of total MR genes`) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  inner_join( # make it the same order as the main matrix
    tibble(
      pheno_label = colnames(msmr_mat),
      order = 1:length(colnames(msmr_mat))
    )
  ) %>%
  arrange(order) %>%
  select(-order) %>%
  column_to_rownames("pheno_label") %>%
  as.matrix()

genes_tested_per_trait_magma$pheno_label == colnames(msmr_mat)

# calculaate the humber of unique genes that are SMR-sig (in at least one cell type) for each trait for the top annotation
n_smr_anno <- df_msmr %>% # this as all sig smr genes
  select(pheno_label, probeID) %>% # get total number of distinct SMRsig genes
  distinct() %>%
  count(pheno_label) %>%
  rename("Number of MR genes" = n) %>%
  inner_join( # make it the same order as the main matrix
    tibble(
      pheno_label = colnames(msmr_mat),
      order = 1:length(colnames(msmr_mat))
    )
  ) %>%
  arrange(order) %>%
  select(-order)

n_smr_anno$pheno_label == colnames(msmr_mat)

# colour function for each cell
col_pal <- paletteer::paletteer_d("MoMAColors::Exter")

col_fun <- colorRamp2(seq(quantile(msmr_mat, 0), quantile(msmr_mat, 1), len = length(col_pal)), col_pal)

discrete_colors <- pal_npg()(5)

# column (top) annotation
top_annotation <- HeatmapAnnotation(
  n_MAGMA = anno_barplot(
    genes_tested_per_trait_magma$`Number of MAGMA significant genes`,
    height = unit(25, "mm"),
    gp = gpar(fill = discrete_colors[4]),
    axis_param = list(side = "right")
  ),
  n_SMR = anno_barplot(
    n_smr_anno$`Number of MR genes`,
    height = unit(25, "mm"),
    gp = gpar(fill = discrete_colors[5]),
    axis_param = list(side = "right")
  ),
  cell_type_specificity = anno_barplot(
    smr_celltype_specificity_mat,
    height = unit(25, "mm"),
    axis_param = list(side = "right"),
    gp = gpar(fill = setNames(pal_viridis(direction = -1)(6), c("1", "2-5", "6-10", "11-15", "16-20", "21-28")))
  ),
  gap = unit(3, "mm"),
  annotation_label = c(
    n_MAGMA = expression(Number ~ of ~ MAGMA ~ eGenes ~ (LFDR[MAGMA] ~ "<" ~ 0.05)),
    n_SMR = "Number of MR genes",
    cell_type_specificity = str_wrap("Proportion of total MR genes", width = 20)
  ),
  annotation_name_rot = 0,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontfamily = "Helvetica")
)

left_annotation <- rowAnnotation(
  NumberOfCells = anno_barplot(log10(df_cell_counts$n_cells),
                               gp = gpar(fill = df_cell_counts$color),
                               height = unit(20, "mm"),
                               axis_param = list(side = "top", direction = "reverse", labels_rot = 90)
  ),
  NumberOfeGenes = anno_barplot(egenes_mat,
                                gp = gpar(fill = c("MR + eGene" = discrete_colors[2], "eGene" = discrete_colors[3])),
                                height = unit(20, "mm"),
                                axis_param = list(side = "top", direction = "reverse", labels_rot = 90)
  ),
  gp = gpar(col = "#404040"),
  gap = unit(3, "mm"),
  annotation_label = c(expression(log[10] ~ (Number ~ of ~ cells)), "Number of eGenes"),
  annotation_name_rot = 45,
  annotation_name_gp = gpar(fontfamily = "Helvetica")
)

# legends for annotations
cell_type_specificity_legend <- Legend(
  at = c("1", "2-5", "6-10", "11-15", "16-20", "21-28"),
  title = str_wrap("Number of cell types with shared MR gene", width = 20),
  type = "grid",
  title_gp = gpar(fontsize = 10, fontface = "plain", fontfamily = "Helvetica"),
  labels_gp = gpar(fontfamily = "Helvetica"),
  legend_gp = gpar(fill = setNames(pal_viridis(direction = -1)(6), c("1", "2-5", "6-10", "11-15", "16-20", "21-28"))),
  direction = "horizontal"
)
number_egenes_legend <- Legend(
  at = c("MR eGene", "eGene"),
  title = str_wrap("Number of eGenes in cell type", width = 20),
  type = "grid",
  title_gp = gpar(fontsize = 10, fontface = "plain", fontfamily = "Helvetica"),
  labels_gp = gpar(fontfamily = "Helvetica"),
  legend_gp = gpar(fill = c("MR eGene" = discrete_colors[2], "eGene" = discrete_colors[3])),
  direction = "horizontal"
)

# colour by decile of the MR gene count rather than the raw count

# Calculate decile breaks from ALL values in the matrix
msmr_vec <- as.vector(msmr_mat)
global_breaks <- quantile(msmr_vec, probs = seq(0, 1, 0.1), na.rm = TRUE) %>%
  round()

# Create labels from the break intervals
break_labels <- paste0(
  c(global_breaks[1], global_breaks[2:(length(global_breaks)-1)] + 1),
  " - ",
  global_breaks[-1]
)

# Apply same breaks to all columns with interval labels
deciles_vec <- quantile(msmr_mat, probs = seq(0, 1, by = 0.1))

decile_matrix <- matrix(
  as.numeric(as.character(cut(msmr_vec, breaks = global_breaks, labels = c(0:9), include.lowest = TRUE))),
  nrow = nrow(msmr_mat),
  ncol = ncol(msmr_mat),
  dimnames = list(rownames(msmr_mat), colnames(msmr_mat))
)

binned_color_pal <- setNames(paletteer::paletteer_d("MoMAColors::Exter"), 0:9)

decile_legend <- Legend(
  at = break_labels,
  title = expression(Number ~ of ~ MR ~ genes ~ (LFDR[MR] ~ "<" ~ 0.05 * "," ~ p[HEIDI] ~ ">=" ~ 0.05)),
  type = "grid",
  title_gp = gpar(fontsize = 10, fontface = "plain", fontfamily = "Helvetica"),
  labels_gp = gpar(fontfamily = "Helvetica"),
  legend_gp = gpar(fill = setNames(binned_color_pal, break_labels)),
  direction = "horizontal"
)

biotrait_heatmap_horizontal <- Heatmap(
  decile_matrix,
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = genes_tested_per_trait_magma$pheno_cat,
  row_split = df_cell_counts$major_cell_type,
  row_title_rot = 0, border = TRUE,
  column_names_rot = 45,
  column_title_rot = 45,
  col = binned_color_pal,
  rect_gp = gpar(col = "white", lwd = 1),
  border_gp = gpar(col = "#161734"),
  row_title_gp = gpar(fontfamily = "Helvetica"),
  column_title_gp = gpar(fontfamily = "Helvetica"),
  row_names_gp = gpar(fontfamily = "Helvetica"),
  column_names_gp = gpar(fontfamily = "Helvetica"),
  left_annotation = left_annotation,
  top_annotation = top_annotation,
  width = ncol(msmr_mat) * unit(5, "mm"),
  height = nrow(msmr_mat) * unit(5, "mm"),
  show_heatmap_legend = FALSE
)

agg_png("figures/v3/supp/mr_heatmap_decile_biotraits.png",
    width = 15 * 2.54, height = 13 * 2.54, units = "cm", res = 300
)
draw(biotrait_heatmap_horizontal,
     annotation_legend_list = list(decile_legend, cell_type_specificity_legend, number_egenes_legend),
     heatmap_legend_side = NULL,
     annotation_legend_side = "right"
)
dev.off()

pdf("figures/v3/supp/mr_heatmap_decile_biotraits.pdf", width = 15, height = 13)
draw(biotrait_heatmap_horizontal,
     annotation_legend_list = list(decile_legend, cell_type_specificity_legend, number_egenes_legend),
     heatmap_legend_side = NULL,
     annotation_legend_side = "right"
)
dev.off()
