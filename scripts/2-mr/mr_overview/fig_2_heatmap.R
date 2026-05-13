# Author: Blake Bowen 
source("scripts/preprocess_strict.R")

suppressPackageStartupMessages({
    library(ragg)
    library(scales)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(glue)
    library(cowplot)
    library(paletteer)
})


# [] create separate similarity heatmap for the biological traits
# [x] plot the two heatmaps together 

# ----------------------------------------------------------------------------------------------------------------------------------
# Plot Parameters 
# ----------------------------------------------------------------------------------------------------------------------------------

PLOT_SIZE <- 0.825  # scale factor: multiply to resize entire plot proportionally
print(glue("Plot size scale factor: {PLOT_SIZE}"))

cell_size_cm <- unit(0.2 * PLOT_SIZE, "cm")
dend_height <- unit(1 * PLOT_SIZE, "cm")
ht_opt$DENDROGRAM_PADDING = unit(0.4 * PLOT_SIZE, "cm") # increase padding to fit in long labels

# ----------------------------------------------------------------------------------------------------------------------------------
# Functions 
# ----------------------------------------------------------------------------------------------------------------------------------

#' calculate jaccard index given two sets of genes 
calculate_jaccard <- function(set_a, set_b) {

    if (is.null(set_a) || is.null(set_b)) {
        return(0)
    }
    union_ab <- length(union(set_a, set_b))
    if (union_ab == 0) {
        return(0)
    }
    jaccard_index <- length(intersect(set_a, set_b)) / union_ab
    # print(jaccard_index)
    return(jaccard_index)
}

#' Simpson index / overlap coefficient quantifies the overlap relative to the smaller set
# ' Simpson index = 1 if one set is a complete subset of another i.e 
#' it penalises the
calculate_simpson_overlap <- function(set_a, set_b) {

    if (is.null(set_a) || is.null(set_b)) {
        return(0)
    }
    # calculate the denominator for the Simpson index (length of the smaller set)
    min_ab <- min(c(length(set_a), length(set_b)))
    if (min_ab == 0) {
        return(0)
    }
    simpson_index <- length(intersect(set_a, set_b)) / min_ab # numerator is the length of the intersect of each set 
    # print(jaccard_index)
    return(simpson_index)
}

# Fisher's exact test
fishers_test_gene_sets <- function(set_a, set_b, set_gene_universe){
    if (is.null(set_a) || is.null(set_b)) {
        return(NA)
    }   
    in_both <- length(intersect(set_a, set_b))

    if (in_both == 0) {
        return(1.0)
    }

    a_only <- length(setdiff(set_a, set_b))
    b_only <- length(setdiff(set_b, set_a))
    in_neither <- length(unique(set_gene_universe)) - (in_both + a_only + b_only)  # this should be all the genes that are not in set a or set b 
    
    contingency_table <- matrix(c(in_both, a_only, b_only, in_neither), nrow=2)
    fisher_test_result <- fisher.test(contingency_table, alternative="greater")
    return(fisher_test_result$p.value)
}

# helper function to convert the tidy data into matrix format for plotting
prepare_ct_matrices <- function(data, value_col, celltype=NULL, is_pvalue = FALSE) {

    # data = plot_data
    # value_col = "jaccard_index"
    # is_pvalue = FALSE

    # 1. filter down to the cell type of interest 
    if (!is.null(celltype)) {
        print("Filtering to keep {celltype}")
        data <- data %>% filter(biosample == {{celltype}})
    }
    
    # 2. Get the lower triangle already calculated 
    df_lower <- data %>% select(p1, p2, value = {{value_col}})
    
    # 3. Create the upper triangle (swap p1 and p2)
    df_upper <- df_lower %>% rename(p1 = p2, p2 = p1)
    
    # 4. Create the diagonal (Self vs Self)
    unique_traits <- unique(c(data$p1, data$p2))
    df_diag <- tibble(
        p1 = unique_traits,
        p2 = unique_traits,
        value = if(is_pvalue) 1 else 1 # P-value 1 (insig) or Jaccard 1 (perfect match)
    )
    
    # 4. Bind and Pivot
    out_matrix <- bind_rows(df_lower, df_upper, df_diag) %>%
        # mutate(row = row_number()) %>%
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

# this plot data is for each cell type specifically
plot_data_trait_level  <- crossing( 
    # Create all pairs of phenotypes within this biosample
        p1 = as.character(phenotypes),
        p2 = as.character(phenotypes)
    ) %>%
    filter(p1 < p2) %>% # Remove duplicates (A vs B is same as B vs A) and self-comparisons so each combo is tested once 
    # Join the gene lists back in
    left_join(nested_genes_trait_level, by = c("p1" = "pheno_label")) %>%
    rename(genes_1 = genes) %>%
    left_join(nested_genes_trait_level, by = c("p2" = "pheno_label")) %>%
    rename(genes_2 = genes) %>% 
    mutate(jaccard_index = map2_dbl(genes_1, genes_2, calculate_jaccard)) %>%  #%>%   # Calculate Jaccard
    mutate(overlap_coefficient = map2_dbl(genes_1, genes_2, calculate_simpson_overlap)) %>% 
    mutate(fishers_test_pvalue = map2_dbl(genes_1, genes_2, fishers_test_gene_sets, set_gene_universe = gene_universe)) %>% 
    # group_by(biosample) %>%
    mutate(
        fishers_test_pvalue_adjusted = p.adjust(fishers_test_pvalue, method = "BH")
    ) %>%
    ungroup() %>%
    mutate()

# add the spearman correlation results (calculated in fig_2_heatmap_rrho_calculation.R)
spearman_corrs <- fread("results/rrho/spearman_corr_all_trait_combos_strictmr.tsv")

plot_data_trait_level <- plot_data_trait_level %>% left_join(spearman_corrs, by = c("p1", "p2"))
plot_data_trait_level %>% filter(is.na(spearman_corr))

# add the RRHO p-values (calculated in fig_2_heatmap_rrho_calculation.R)
rrho_pval <- fread("results/rrho/rrho_pvalues_all_trait_combos_two.sided.tsv")

plot_data_trait_level <- plot_data_trait_level %>%
    left_join(rrho_pval, by = c("p1", "p2")) %>% 
    # reverse -log10 transformation to get the BY p-values
    mutate(max_p_value_BY = 10^(-max_p_value_BY))



# ----------------------------------------------------------------------------------------------------------------------------------
# generate matrices to plot in the main heatmap 
# ----------------------------------------------------------------------------------------------------------------------------------

# make the full matrices
mat_spearman <- prepare_ct_matrices(plot_data_trait_level, spearman_corr, is_pvalue = FALSE)
mat_pvalue <- prepare_ct_matrices(plot_data_trait_level, max_p_value_BY, is_pvalue = TRUE)
mat_pvalue <- mat_pvalue[rownames(mat_spearman), colnames(mat_spearman)]

# mat_jaccard <- prepare_ct_matrices(data = plot_data_trait_level, value_col = jaccard_index, is_pvalue = FALSE)
# mat_overlap_coef <- prepare_ct_matrices(data = plot_data_trait_level, value_col = overlap_coefficient, is_pvalue = FALSE)
# mat_pvalue  <- prepare_ct_matrices(plot_data_trait_level, fishers_test_pvalue_adjusted, is_pvalue = TRUE)
# mat_pvalue  <- mat_pvalue[rownames(mat_jaccard), colnames(mat_jaccard)]

# Get the trait metadata to include as annotation bars
trait_meta <- fread("resources/metadata/trait_metadata_n.tsv") %>% 
    filter(include==TRUE)
# source("scripts/mr_overview/temp.R") # run this way until the file above is fixed 
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
# disease_mat_pvalue <- mat_pvalue[disease_traits,disease_traits]

bio_mat_spearman <- mat_spearman[bio_traits,bio_traits]
# bio_mat_pvalue <- mat_pvalue[bio_traits,bio_traits]

# disease_mat_overlap_coef <- mat_overlap_coef[disease_traits,disease_traits]
# bio_mat_overlap_coef <- mat_overlap_coef[bio_traits,bio_traits]
# TODO: 
# [] There are a lot of traits in the metadata that don't appear to be in the results... will remove for now but need to investigate further why these are missing / if they are suppposed to be in there .. 

# look at the n intersecting and union sizes 
# not actually used for plotting right now
# plot_data_bio <- plot_data_trait_level %>%
# filter(
#     p1 %in% bio_traits,
#     p2 %in% bio_traits
# )
# bio_traits_data <- plot_data_bio %>%
#     mutate(
#         intersecting = map2(genes_1, genes_2, intersect),
#         union = map2(genes_1, genes_2, union),
#     ) %>%
#     rowwise() %>% 
#     mutate(
#         n_intersecting = length(intersecting),
#         n_union =  length(union),
#         n_genes_1 = length(genes_1),
#         n_genes_2 = length(genes_2),       
#     ) %>% 
#     select(
#         -intersecting, -union
#     )

# bio_traits_data %>% filter(p1 == "eBmd" | p2 == "eBmd")  %>% select(p1, p2, n_genes_1, n_genes_2, n_intersecting, n_union, jaccard_index, overlap_coefficient) %>% write_tsv("temp/embd_jaccard_summary.tsv")


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
    annotation_name_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica")
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
    annotation_name_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica")
)

# ----------------------------------------------------------------------------------------------------------------------------------
# Disease traits  -  Plot the heatmap for all cell types 
# ----------------------------------------------------------------------------------------------------------------------------------

dend_height <- unit(1 * PLOT_SIZE, "cm") # set dendrogram height

legend_disease  <- packLegend(
        Legend(
            col_fun = col_fun,
            title = expression("Spearman's" ~ rho),
            title_gp = gpar(fontsize = 8 * PLOT_SIZE, fontfamily = "Helvetica"),
            labels_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),
            grid_height = unit(0.3 * PLOT_SIZE, "cm"),
            grid_width = unit(0.3 * PLOT_SIZE, "cm"),
            legend_height = unit(length(disease_category_cols) * 0.3 * PLOT_SIZE, "cm"),
        ),
        Legend(
             labels = names(disease_category_cols),
             legend_gp = gpar(fill = disease_category_cols),
             title = "Trait category",
             title_gp = gpar(fontsize = 8 * PLOT_SIZE, fontfamily = "Helvetica"),
             labels_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),
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
            
            # p1 <- rownames(disease_mat_jaccard[i, j])
            # p2 <- colnames(disease_mat_jaccard[i, j])

            # plot rectangles if the phenotypes are different (ie. not the diagonal) 
            if (i != j){
                # or if i == j?
                grid.rect(x, y, w, h, 
                    gp = gpar(fill = fill, col = "grey", lwd = 0.5))
            } else {
                grid.rect(x, y, w, h,
                    gp = gpar(col = "grey", fill = "transparent", lwd = 0.5))
            }

            # plot significance 
        #     p_val <- mat_pvalue[i, j]
        #     if(p_val < 0.05 ){
        #          grid.text("*", x, y, gp = gpar(fontsize = 15))
        #     }
		}
    },

    # Annotations:
    bottom_annotation = disease_bottom_ha,

    # Labels
    column_title = " ",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 0),
    row_names_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),
    row_names_side = "left",
    row_labels = disease_trait_meta$name,    # Use nice names for rows
    column_labels = disease_trait_meta$name,
    column_names_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),

    # Dendrograms
    column_dend_side = "bottom",
    row_dend_width = dend_height,
    column_dend_height = dend_height,
    # show_row_dend = FALSE,
    show_column_dend = FALSE,

    # Legend
)

dend_height <- unit(1 * PLOT_SIZE, "cm") # set dendrogram height

legend_bio <- packLegend(
        Legend(
            col_fun = col_fun,
            title = expression("Spearman's" ~ rho),
            title_gp = gpar(fontsize = 8 * PLOT_SIZE, fontfamily = "Helvetica"),
            labels_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),
            grid_height = unit(0.3 * PLOT_SIZE, "cm"),
            grid_width = unit(0.3 * PLOT_SIZE, "cm"),
            legend_height = unit(length(bio_category_cols) * 0.3 * PLOT_SIZE, "cm"),
        ),
        Legend(
            title = "Trait category",
            labels = names(bio_category_cols),
            legend_gp = gpar(fill = bio_category_cols),
            title_gp = gpar(fontsize = 8 * PLOT_SIZE, fontfamily = "Helvetica"),
            labels_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),
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
            
            # p1 <- rownames(bio_mat_jaccard[i, j])
            # p2 <- colnames(bio_mat_jaccard[i, j])

            # plot rectangles if the phenotypes are different (ie. not the diagonal) 
            if (i != j){
                # or if i == j?
                grid.rect(x, y, w, h, 
                    gp = gpar(fill = fill, col = "grey", lwd = 0.5))
            } else {
                grid.rect(x, y, w, h,
                    gp = gpar(col = "grey", fill = "transparent", lwd = 0.5))
            }

            # plot significance 
            # p_val <- mat_pvalue[i, j]
            # if(p_val < 0.05 ){
            #      grid.text("*", x, y, gp = gpar(fontsize = 15))
            # }
		}
    },

    # Annotations:
    bottom_annotation = bio_bottom_ha,

    # Labels
    column_title = " ",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 0),
    row_names_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),
    row_names_side = "left",
    row_labels = bio_trait_meta$name,    # Use nice names for rows
    column_labels = bio_trait_meta$name,
    column_names_gp = gpar(fontsize = 7 * PLOT_SIZE, fontfamily = "Helvetica"),

    # Dendrograms
    column_dend_side = "bottom",
    row_dend_width = dend_height,
    column_dend_height = dend_height,
    # show_row_dend = FALSE,
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

pdf(glue("fig_pub/figure_2_heatmap_all_cell_types_disease_only.pdf"), width = pdf_width_in, height = pdf_height_in)

draw(disease_heatmap, show_heatmap_legend = FALSE)
decorate_heatmap_body("Spearman", {
    draw(legend_disease, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
})
dev.off()

disease_heatmap %>% saveRDS("fig_pub/figure_2_heatmap_disease_only.rds")
legend_disease %>% saveRDS("fig_pub/figure_2_heatmap_legend_disease.rds")


# ----------------------------------------------------------------------------------------------------------------------------------
# Draw separate disease  and biologic trait heatmaps
# ----------------------------------------------------------------------------------------------------------------------------------


# heatmap_list <- bio_heatmap + disease_heatmap

pdf_width_in <- as.numeric((ncol(disease_mat_spearman + disease_mat_spearman) * cell_size_cm) + unit(35 * PLOT_SIZE, "cm"))/2.54
pdf_height_in <- as.numeric((nrow(disease_mat_spearman) * cell_size_cm) + unit(10 * PLOT_SIZE, "cm"))/2.54

ht_opt$DENDROGRAM_PADDING = unit(1.05 * PLOT_SIZE, "cm") # increase padding to fit in long labels
p1 = grid.grabExpr(draw(disease_heatmap))
ht_opt$DENDROGRAM_PADDING = unit(2.1 * PLOT_SIZE, "cm") # increase padding to fit in long labels
p2 = grid.grabExpr(draw(bio_heatmap))

pdf(glue("fig_pub/figure_2_heatmap_all_cell_types_disease_and_biological_separate.pdf"), width = pdf_width_in, height = pdf_height_in)

cowplot::plot_grid(p1, p2)

dev.off()





####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# EXPERIMENTAL: SIMPSON INDEX HEATMAPS 

# col_fun_simpson <- colorRamp2(
#     breaks = seq(0, max(mat_overlap_coef), length.out = 9), 
#     colors = brewer.pal(9, "YlOrRd")
# )


# dend_height <- unit(3, "cm") # set dendrogram height 

# disease_heatmap_simpson <- Heatmap(
#     disease_mat_overlap_coef, 
#     name = "Jaccard",

#     # colour  
#     # col = disease_col_fun,
#     col = col_fun_simpson,
#     # cell dimensions 
#     width  = cell_size_cm * ncol(disease_mat_jaccard),
#     height = cell_size_cm * nrow(disease_mat_jaccard),

#     # Draw Lower Triangle Only
#     rect_gp = gpar(type = "none"), # wipe default graphic params 
#     cell_fun = function(j, i, x, y, w, h, fill) {

#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            
#             # p1 <- rownames(disease_mat_jaccard[i, j])
#             # p2 <- colnames(disease_mat_jaccard[i, j])

#             # plot rectangles if the phenotypes are different (ie. not the diagonal) 
#             if (i != j){
#                 # or if i == j?
#                 grid.rect(x, y, w, h, 
#                     gp = gpar(fill = fill, col = fill))
#             } else {
#                 grid.rect(x, y, w, h,
#                     gp = gpar(col = "grey", fill = "transparent", lwd = 0.5))
#             }

#             # plot fisher's test significance 
#             # p_val <- mat_pvalue[i, j]
#             # if(p_val < 0.05 ){
#             #      grid.text("*", x, y, gp = gpar(fontsize = 15))
#             # }
# 		}
#     },

#     # Annotations:
#     left_annotation = disease_row_ha,
        
#     # Labels
#     column_title = paste("Overlap Coefficient: disease traits"),
#     row_names_gp = gpar(fontsize = 14),
#     row_names_side = "left",
#     row_labels = disease_trait_meta$name,    # Use nice names for rows
#     column_labels = disease_trait_meta$name,
#     column_names_gp = gpar(fontsize = 14),

#     # Dendrograms
#     column_dend_side = "bottom",
#     row_dend_width = dend_height,
#     column_dend_height = dend_height,
#     # show_row_dend = FALSE,
#     # show_column_dend = FALSE
# )

# dend_height <- unit(5, "cm") # set dendrogram height 

# bio_heatmap_simpson <- Heatmap(
#     bio_mat_overlap_coef, 
#     name = "Overlap Coefficient",

#     # colour  
#     col = col_fun_simpson,
#     # cell dimensions 
#     width  = cell_size_cm * ncol(bio_mat_jaccard),
#     height = cell_size_cm * nrow(bio_mat_jaccard),

#     # Draw Lower Triangle Only
#     rect_gp = gpar(type = "none"), # wipe default graphic params 
#     cell_fun = function(j, i, x, y, w, h, fill) {

#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            
#             # p1 <- rownames(bio_mat_jaccard[i, j])
#             # p2 <- colnames(bio_mat_jaccard[i, j])

#             # plot rectangles if the phenotypes are different (ie. not the diagonal) 
#             if (i != j){
#                 # or if i == j?
#                 grid.rect(x, y, w, h, 
#                     gp = gpar(fill = fill, col = fill))
#             } else {
#                 grid.rect(x, y, w, h,
#                     gp = gpar(col = "grey", fill = "transparent", lwd = 0.5))
#             }

#             # plot fisher's test significance 
#             # p_val <- mat_pvalue[i, j]
#             # if(p_val < 0.05 ){
#             #      grid.text("*", x, y, gp = gpar(fontsize = 15))
#             # }
# 		}
#     },

#     # Annotations:
#     left_annotation = bio_row_ha,
        
#     # Labels
#     column_title = paste("Overlap Coefficient: biological traits"),
#     row_names_gp = gpar(fontsize = 14),
#     row_names_side = "left",
#     row_labels = bio_trait_meta$name,    # Use nice names for rows
#     column_labels = bio_trait_meta$name,
#     column_names_gp = gpar(fontsize = 14),

#     # Dendrograms
#     column_dend_side = "bottom",
#     row_dend_width = dend_height,
#     column_dend_height = dend_height,
#     # show_row_dend = FALSE,
#     # show_column_dend = FALSE
# )


# # ----------------------------------------------------------------------------------------------------------------------------------
# # SIMPSON Draw separate disease  and biologic trait heatmaps
# # ----------------------------------------------------------------------------------------------------------------------------------

# # heatmap_list <- bio_heatmap + disease_heatmap

# pdf_width_in <- as.numeric((ncol(disease_mat_overlap_coef + disease_mat_overlap_coef) * cell_size_cm) + unit(70, "cm"))/2.54
# pdf_height_in <- as.numeric((nrow(disease_mat_overlap_coef) * cell_size_cm) + unit(20, "cm"))/2.54

# ht_opt$DENDROGRAM_PADDING = unit(3, "cm") # increase padding to fit in long labels 
# p1 = grid.grabExpr(draw(disease_heatmap_simpson))
# ht_opt$DENDROGRAM_PADDING = unit(6, "cm") # increase padding to fit in long labels 
# p2 = grid.grabExpr(draw(bio_heatmap_simpson))

# pdf(glue("fig_pub/figure_2_heatmap_all_cell_types_disease_and_biological_separate_simpson_index.pdf"), width = pdf_width_in, height = pdf_height_in)

# cowplot::plot_grid(p1, p2)

# dev.off()


# ####################################################################################################################################
# ####################################################################################################################################
# ####################################################################################################################################

# # ----------------------------------------------------------------------------------------------------------------------------------
# # get gene sets for each trait x cell type 
# # ----------------------------------------------------------------------------------------------------------------------------------

# # cell type specific df 
# nested_genes <- df_msmr_strict %>%
#   select(biosample, pheno_label, probeID) %>%
#   group_by(biosample, pheno_label) %>%
#   summarise(genes = list(unique(probeID)), .groups = "drop") %>% 
#   ungroup()


# # ----------------------------------------------------------------------------------------------------------------------------------
# # individual cell types heatmap plot 
# # ----------------------------------------------------------------------------------------------------------------------------------


# # this plot data is for each cell type specifically
# plot_data <- nested_genes %>%
#     group_by(biosample) %>%
#     # Create all pairs of pheno_labels within this biosample
#     reframe(crossing(
#         p1 = pheno_label, 
#         p2 = pheno_label
#     )) %>%
#     filter(p1 < p2) %>% # Remove duplicates (A vs B is same as B vs A) and self-comparisons so each combo is tested once 
#     # Join the gene lists back in
#     left_join(nested_genes, by = c("biosample", "p1" = "pheno_label")) %>%
#     rename(genes_1 = genes) %>%
#     left_join(nested_genes, by = c("biosample", "p2" = "pheno_label")) %>%
#     rename(genes_2 = genes) %>% 
#     mutate(jaccard_index = map2_dbl(genes_1, genes_2, calculate_jaccard)) %>%  #%>%   # Calculate Jaccard
#     mutate(fishers_test_pvalue = map2_dbl(genes_1, genes_2, fishers_test_gene_sets, set_gene_universe = gene_universe)) %>% 
#     group_by(biosample) %>%
#     mutate(
#         # OPTION B: Stratified Adjustment (Per Cell Type)
#         # Only use this if each cell type is a totally separate hypothesis/paper section
#         fishers_test_pvalue_adjusted = p.adjust(fishers_test_pvalue, method = "BH")
#     ) %>%
#     ungroup() %>% 
#     mutate()

# # ----------------------------------------------------------------------------------------------------------------------------------
# # generate matrices to plot in the main heatmap 
# # ----------------------------------------------------------------------------------------------------------------------------------

# celltype_plot <- "B_intermediate"
# mat_jaccard <- prepare_ct_matrices(plot_data, jaccard_index, celltype_plot, is_pvalue = FALSE)
# mat_pvalue  <- prepare_ct_matrices(plot_data, fishers_test_pvalue_adjusted, celltype_plot, is_pvalue = TRUE)
# mat_pvalue  <- mat_pvalue[rownames(mat_jaccard), colnames(mat_jaccard)]

# # Get the trait metadata to include as annotation bars
# # trait_meta <- fread("resources/misc/trait_metadata_n.tsv")
# source("scripts/mr_overview/temp.R") # run this way until the file above is fixed 
# trait_meta <- trait_meta %>% as.data.frame() %>% column_to_rownames("trait_id")
# trait_meta_ordered <- trait_meta[rownames(mat_jaccard), ]

# col_fun <- colorRamp2(
#     breaks = seq(0, max(mat_jaccard), length.out = 9), 
#     colors = brewer.pal(9, "YlOrRd")
# )

# # col_fun <- colorRamp2(
# #     breaks = seq(0, max(mat_jaccard), length.out = 11), 
# #     colors = rev(brewer.pal(11, "RdYlBu"))
# # )


# # ----------------------------------------------------------------------------------------------------------------------------------

# ht <- Heatmap(
#     mat_jaccard, 
#     name = "Jaccard",

#     # colour  
#     col = col_fun,

#     width  = cell_size_cm * ncol(mat_jaccard),
#     height = cell_size_cm * nrow(mat_jaccard),

#     # draw asterisks 
#     # cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#     #     grid.text(mat_pvalue[i, j], x, y)
#     # },

#     rect_gp = gpar(type = "none"),

#     # Draw Lower Triangle Only
#     cell_fun = function(j, i, x, y, w, h, fill) {

#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             # plot rectangles 
# 			grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#             # plot fisher's test significance 
#             p_val <- mat_pvalue[i, j]
#             if(p_val < 0.05 ){
#                  grid.text("*", x, y, gp = gpar(fontsize = 15))
#             }
# 		}
#     },

#     # Annotations:
#     left_annotation = row_ha,
        
#     # Labels
#     column_title = paste("Jaccard Overlap:", celltype_plot),
#     row_names_gp = gpar(fontsize = 14),
#     row_names_side = "left",
#     row_labels = trait_meta_ordered$name,    # Use nice names for rows
#     column_labels = trait_meta_ordered$name,
#     column_names_gp = gpar(fontsize = 14),

#     # dendrograms
#     column_dend_side = "bottom",
#     row_dend_width = dend_height,
#     column_dend_height = dend_height,
#     # show_row_dend = FALSE,
#     # show_column_dend = FALSE
# )


# # pdf height and width 
# pdf_width_in <- as.numeric((ncol(mat_jaccard) * cell_size_cm) + unit(30, "cm"))/2.54
# pdf_height_in <- as.numeric((nrow(mat_jaccard) * cell_size_cm) + unit(20, "cm"))/2.54

# pdf(glue("fig_pub/figure_2_heatmap_{celltype_plot}_with_p.pdf"), width = pdf_width_in, height = pdf_height_in)

# draw(ht)

# dev.off()


