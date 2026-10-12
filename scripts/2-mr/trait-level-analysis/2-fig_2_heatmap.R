# Author: Blake Bowen
source("scripts/0-preprocess/preprocess_results.R")

suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(glue)
    library(cowplot)
    library(paletteer)
})

# fixed sort order so trait pairs (p1 < p2) match the saved correlation tables regardless of IDE/terminal locale
Sys.setlocale("LC_COLLATE", "C")

dir.create("figures/mr_overview", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------------------------------------------------------------
# Plot Parameters
# ----------------------------------------------------------------------------------------------------------------------------------

PLOT_SIZE <- 0.7  # scale factor: multiply to resize entire plot proportionally
print(glue("Plot size scale factor: {PLOT_SIZE}"))

cell_size_cm <- unit(0.2 * PLOT_SIZE, "cm")
dend_height <- unit(1 * PLOT_SIZE, "cm")
ht_opt$DENDROGRAM_PADDING = unit(0.4 * PLOT_SIZE, "cm") # increase padding to fit in long labels

# ----------------------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------------------

# helper function to convert the tidy pairwise data into a symmetric matrix (1 on the diagonal) for plotting
prepare_ct_matrices <- function(data, value_col) {

    # 1. Get the lower triangle already calculated
    df_lower <- data %>% select(p1, p2, value = {{value_col}})

    # 2. Create the upper triangle (swap p1 and p2)
    df_upper <- df_lower %>% rename(p1 = p2, p2 = p1)

    # 3. Create the diagonal (Self vs Self)
    unique_traits <- unique(c(data$p1, data$p2))
    df_diag <- tibble(
        p1 = unique_traits,
        p2 = unique_traits,
        value = 1
    )

    # 4. Bind and Pivot
    out_matrix <- bind_rows(df_lower, df_upper, df_diag) %>%
        pivot_wider(names_from = p2, values_from = value, values_fill = NA) %>%
        column_to_rownames("p1") %>%
        as.matrix()

    out_matrix <- out_matrix[,rownames(out_matrix)]

    return(out_matrix)
}

# ----------------------------------------------------------------------------------------------------------------------------------
# get significant gene sets for each trait
# ----------------------------------------------------------------------------------------------------------------------------------

df_msmr <- df_msmr_tenk10k %>% filter(mr)

# all cell types combined
# take the set of siginficant genes for each trait across all cell types
nested_genes_trait_level <-  df_msmr %>%
    select(pheno_label, probeID) %>%
    group_by(pheno_label) %>%
    summarise(genes = list(unique(probeID)), .groups = "drop") %>%
    ungroup() %>%
    rowwise() %>%
    mutate(n_genes = length(genes)) %>%
    ungroup() %>%
    filter(n_genes >= 10) %>%
    select(-n_genes)

phenotypes <- unique(nested_genes_trait_level$pheno_label)

plot_data_trait_level  <- crossing(
        p1 = as.character(phenotypes),
        p2 = as.character(phenotypes)
    ) %>%
    filter(p1 < p2) # Remove duplicates (A vs B is same as B vs A) and self-comparisons so each combo is tested once

# add the spearman correlation results
spearman_corrs <- fread("results/rrho/spearman_corr_all_trait_combos_strictmr.tsv")

plot_data_trait_level <- plot_data_trait_level %>% left_join(spearman_corrs, by = c("p1", "p2"))
plot_data_trait_level %>% filter(is.na(spearman_corr))

# ----------------------------------------------------------------------------------------------------------------------------------
# generate matrices to plot in the main heatmap
# ----------------------------------------------------------------------------------------------------------------------------------

mat_spearman <- prepare_ct_matrices(plot_data_trait_level, spearman_corr)

# Get the trait metadata to include as annotation bars
trait_meta <- fread("metadata/trait.tsv") %>%
    filter(include==TRUE)
trait_meta <- trait_meta %>% as.data.frame() %>% column_to_rownames("label")
trait_meta_ordered <- trait_meta[rownames(mat_spearman), ]

# split biological vs disease categories into separate heatmaps
bio_trait_meta <- trait_meta_ordered %>%
    filter(supercategory == "biological")
disease_trait_meta <- trait_meta_ordered %>%
    filter(supercategory == "disease")

disease_traits <- rownames(disease_trait_meta)
bio_traits <- rownames(bio_trait_meta)

disease_mat_spearman <- mat_spearman[disease_traits,disease_traits]
bio_mat_spearman <- mat_spearman[bio_traits,bio_traits]

# --------------------------------------------------------------------------------------------------------------------------------------
# Color mappings
# --------------------------------------------------------------------------------------------------------------------------------------

# color function for spearman corellations
col_fun <- colorRamp2(
    breaks = seq(-0.7, 0.7, length.out = 9),
    colors = rev(brewer.pal(9, "RdYlBu"))
)

# --------------------------------------------------------------------------------------------------------------------------------------
#  disease traits  - Create the annotations
# ----------------------------------------------------------------------------------------------------------------------------------

trait_cat_col <- paletteer_d("ggthemes::Tableau_10", 10) %>%
  set_names(cat_order)

disease_category_cols <- setNames(trait_cat_col[which(names(trait_cat_col) %in% unique(disease_trait_meta$cat_rev))], names(trait_cat_col)[which(names(trait_cat_col) %in% unique(disease_trait_meta$cat_rev))])
bio_category_cols <- setNames(trait_cat_col[which(names(trait_cat_col) %in% unique(bio_trait_meta$cat_rev))], names(trait_cat_col)[which(names(trait_cat_col) %in% unique(bio_trait_meta$cat_rev))])


disease_bottom_ha <- HeatmapAnnotation(
    Category = disease_trait_meta$cat_rev,
    gp = gpar(col = "grey", lwd = 0.5), # add border to annotation bars
    simple_anno_size = cell_size_cm,
    col = list(
        Category = disease_category_cols,
        labels = names(disease_category_cols)
    ),
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 5, fontfamily = "Helvetica")
)

bio_bottom_ha <- HeatmapAnnotation(
    Category = bio_trait_meta$cat_rev,
    gp = gpar(col = "grey", lwd = 0.5), # add border to annotation bars
    simple_anno_size = cell_size_cm,
    col = list(
        Category = bio_category_cols,
        labels = names(bio_category_cols)
    ),
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 5, fontfamily = "Helvetica")
)

# ----------------------------------------------------------------------------------------------------------------------------------
# Disease traits  -  Plot the heatmap for all cell types
# ----------------------------------------------------------------------------------------------------------------------------------

dend_height <- unit(1 * PLOT_SIZE, "cm") # set dendrogram height

legend_disease  <- packLegend(
        Legend(
            col_fun = col_fun,
            title = expression("Spearman's" ~ rho),
            title_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
            labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
            grid_height = unit(0.3 * PLOT_SIZE, "cm"),
            grid_width = unit(0.3 * PLOT_SIZE, "cm"),
            legend_height = unit(length(disease_category_cols) * 0.3 * PLOT_SIZE, "cm"),
        ),
        Legend(
             labels = names(disease_category_cols),
             legend_gp = gpar(fill = disease_category_cols),
             title = "Trait category",
             title_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
             labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
             grid_height = unit(0.3 * PLOT_SIZE, "cm"),
             grid_width = unit(0.3 * PLOT_SIZE, "cm"),
        ),
        direction = "horizontal"
)

disease_heatmap <- Heatmap(
    disease_mat_spearman,
    name = "Spearman",

    # colour
    col = col_fun,
    # cell dimensions
    width  = cell_size_cm * ncol(disease_mat_spearman),
    height = cell_size_cm * nrow(disease_mat_spearman),

    # Draw Lower Triangle Only
    rect_gp = gpar(type = "none"), # wipe default graphic params
    cell_fun = function(j, i, x, y, w, h, fill) {

        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {

            # plot rectangles if the phenotypes are different (ie. not the diagonal)
            if (i != j){
                grid.rect(x, y, w, h,
                    gp = gpar(fill = fill, col = "grey", lwd = 0.5))
            } else {
                grid.rect(x, y, w, h,
                    gp = gpar(col = "grey", fill = "transparent", lwd = 0.5))
            }
		}
    },

    # Annotations:
    bottom_annotation = disease_bottom_ha,

    # Labels
    column_title = character(0), # a blank " " title still reserves ~5 mm of padding below the column names
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 0),
    row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
    row_names_side = "left",
    row_labels = disease_trait_meta$name,    # Use nice names for rows
    column_labels = disease_trait_meta$name,
    column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),

    # Dendrograms
    column_dend_side = "bottom",
    row_dend_width = dend_height,
    column_dend_height = dend_height,
    show_row_dend = FALSE,
    show_column_dend = FALSE
)

dend_height <- unit(1 * PLOT_SIZE, "cm") # set dendrogram height

legend_bio <- packLegend(
        Legend(
            col_fun = col_fun,
            title = expression("Spearman's" ~ rho),
            title_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
            labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
            grid_height = unit(0.3 * PLOT_SIZE, "cm"),
            grid_width = unit(0.3 * PLOT_SIZE, "cm"),
            legend_height = unit(length(bio_category_cols) * 0.3 * PLOT_SIZE, "cm"),
        ),
        Legend(
            title = "Trait category",
            labels = names(bio_category_cols),
            legend_gp = gpar(fill = bio_category_cols),
            title_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
            labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
            grid_height = unit(0.3 * PLOT_SIZE, "cm"),
            grid_width = unit(0.3 * PLOT_SIZE, "cm"),
        ),
        direction = "horizontal"
)

bio_heatmap <- Heatmap(
    bio_mat_spearman,
    name = "Spearman",

    # colour
    col = col_fun,
    # cell dimensions
    width  = cell_size_cm * ncol(bio_mat_spearman),
    height = cell_size_cm * nrow(bio_mat_spearman),

    # Draw Lower Triangle Only
    rect_gp = gpar(type = "none"), # wipe default graphic params
    cell_fun = function(j, i, x, y, w, h, fill) {

        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {

            # plot rectangles if the phenotypes are different (ie. not the diagonal)
            if (i != j){
                grid.rect(x, y, w, h,
                    gp = gpar(fill = fill, col = "grey", lwd = 0.5))
            } else {
                grid.rect(x, y, w, h,
                    gp = gpar(col = "grey", fill = "transparent", lwd = 0.5))
            }
		}
    },

    # Annotations:
    bottom_annotation = bio_bottom_ha,

    # Labels
    column_title = " ",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 0),
    row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
    row_names_side = "left",
    row_labels = bio_trait_meta$name,    # Use nice names for rows
    column_labels = bio_trait_meta$name,
    column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),

    # Dendrograms
    column_dend_side = "bottom",
    row_dend_width = dend_height,
    column_dend_height = dend_height,
    show_column_dend = FALSE,

    # Legend
    show_heatmap_legend = FALSE
)

# ----------------------------------------------------------------------------------------------------------------------------------
# Draw heatmap - disease only
# ----------------------------------------------------------------------------------------------------------------------------------

# pdf height and width
pdf_width_in <- as.numeric((ncol(disease_mat_spearman) * cell_size_cm) + unit(7 * PLOT_SIZE, "cm"))/2.54
pdf_height_in <- as.numeric((nrow(disease_mat_spearman) * cell_size_cm) + unit(6 * PLOT_SIZE, "cm"))/2.54

pdf(glue("figures/mr_overview/figure_2_heatmap_all_cell_types_disease_only.pdf"), width = pdf_width_in, height = pdf_height_in)

draw(disease_heatmap, show_heatmap_legend = FALSE)
decorate_heatmap_body("Spearman", {
    draw(legend_disease, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
})
dev.off()

disease_heatmap %>% saveRDS("figures/mr_overview/figure_2_heatmap_disease_only.rds")
legend_disease %>% saveRDS("figures/mr_overview/figure_2_heatmap_legend_disease.rds")


# ----------------------------------------------------------------------------------------------------------------------------------
# Draw separate disease  and biologic trait heatmaps
# ----------------------------------------------------------------------------------------------------------------------------------

pdf_width_in <- as.numeric((ncol(disease_mat_spearman + disease_mat_spearman) * cell_size_cm) + unit(35 * PLOT_SIZE, "cm"))/2.54
pdf_height_in <- as.numeric((nrow(disease_mat_spearman) * cell_size_cm) + unit(10 * PLOT_SIZE, "cm"))/2.54

ht_opt$DENDROGRAM_PADDING = unit(1.05 * PLOT_SIZE, "cm") # increase padding to fit in long labels
p1 = grid.grabExpr(draw(disease_heatmap))
ht_opt$DENDROGRAM_PADDING = unit(2.1 * PLOT_SIZE, "cm") # increase padding to fit in long labels
p2 = grid.grabExpr(draw(bio_heatmap))

pdf(glue("figures/mr_overview/figure_2_heatmap_all_cell_types_disease_and_biological_separate.pdf"), width = pdf_width_in, height = pdf_height_in)

cowplot::plot_grid(p1, p2)

dev.off()
