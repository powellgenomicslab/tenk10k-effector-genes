# Author: Blake Bowen 
# mm activate causal-env

library(patchwork)
library(glue)
library(janitor)
library(tidyverse)
library(arrow)
library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)

# source("scripts/mr_overview/rrho_function.R")

# ----------------------------------------------------------------------------------------------------------------------------------
# Plot Parameters 
# ----------------------------------------------------------------------------------------------------------------------------------

cell_size_cm <- unit(0.4, "cm")
dend_height <- unit(3, "cm")
ht_opt$DENDROGRAM_PADDING = unit(3, "cm") # increase padding to fit in long labels 

# ----------------------------------------------------------------------------------------------------------------------------------
# Functions 
# ----------------------------------------------------------------------------------------------------------------------------------

calc_spearman <- function(p1, p2){
    print(glue("Calculating spearman correllation for {p1} vs {p2}..."))

    df1 <- trait_lvl_genes_ranked_by_pSMR_multi %>%
        filter(pheno_label == p1) %>%
        select(probeID, median_rank_normalised_reranked)
    
    df2 <- trait_lvl_genes_ranked_by_pSMR_multi %>%
        filter(pheno_label == p2) %>%
        select(probeID, median_rank_normalised_reranked)

    # Use genes that were tested in both traits
    intersecting_genes <- intersect(df1$probeID, df2$probeID)
    df1 <- df1 %>% filter(probeID %in% intersecting_genes)
    # print(head(df1))

    df2 <- df2 %>% filter(probeID %in% intersecting_genes)
    # print(head(df2))

    print(glue("Spearman correlation on {length(intersecting_genes)} genes that were tested in both traits."))

    # calculate spearman correlation between the two ranked lists
    spearman_corr <- cor(
        df1$median_rank_normalised_reranked,
        df2[match(df1$probeID, df2$probeID)]$median_rank_normalised_reranked,
        method = "spearman"
    )

    tibble(
            p1 = p1,
            p2 = p2,
            spearman_corr = spearman_corr,
            n_intersecting_genes = length(intersecting_genes)
        ) %>% return()

}

# ----------------------------------------------------------------------------------------------------------------------------------
#  Data preprocessing  
# ----------------------------------------------------------------------------------------------------------------------------------

# source("scripts/0-preprocess/preprocess_results.R")
df_msmr_strict <- read_parquet("/g/data/fy54/analysis/tenk10k-causal/results/preprocessed/tenk10k_phase1.v4.parquet.gz") %>%
    setDT()

# start with the full mr result list - not filtered for sig only 
#  use strict msmr results, which are stringently filtered to only include more significant eGenes
pheno_order <- df_msmr_strict[, .N, by = pheno_label][order(-N), pheno_label]
df_msmr_strict[, pheno_label := factor(pheno_label, pheno_order)]
df_msmr_strict <- df_msmr_strict[p_SMR_multi < 0.05,] # filter to keep only nominally significant genes to reduce noise 

# try out some approach of filtering to only GWAS-significant genes?
# df_msmr_strict <- df_msmr_strict[p_GWAS < 5e-8]


# ----------------------------------------------------------------------------------------------------------------------------------
# get significant gene sets for each trait 
# ----------------------------------------------------------------------------------------------------------------------------------

# all cell types combined 
# method for combining gene lists across cell types, to get a single ranked gene list per trait:
# 1) rank genes by signed -log10(p_SMR_multi) for each trait and each cell type, take the median rank across cell types

# TODO: try correllation of the betas instead

df_msmr_strict[, .(n_genes_per_celltype = .N), by = .(pheno_label, biosample)] %>%
    arrange(n_genes_per_celltype) 

trait_lvl_genes_ranked_by_pSMR_multi <- df_msmr_strict[,
        .(  
            probeID,
            p_SMR_multi,
            b_SMR,
            pheno_label,
            biosample,
            signed_log10_p_SMR_multi = (-log10(p_SMR_multi))*sign(b_SMR)
        ),
    ][, 
        .(
            probeID,
            p_SMR_multi,
            b_SMR,
            signed_log10_p_SMR_multi,
            rank = frank(signed_log10_p_SMR_multi, ties.method = "average"),
            normalised_rank = frank(signed_log10_p_SMR_multi, ties.method = "average") / .N # normalised rank between 0 and 1 (normalised by number of genes per cell type)
        ), 
        by = .(pheno_label, biosample) 
    ][,
        .(
            median_rank = median(rank),
            median_rank_normalised = median(normalised_rank),
            # mean_rank = mean(normalised_rank),
            all_signed_log10_p_SMR_multi = list(signed_log10_p_SMR_multi),
            all_ranks = list(normalised_rank),
            median_signed_log10_p_SMR_multi = median(signed_log10_p_SMR_multi),
            all_cell_types = list(biosample),
            n_cell_types = .N
    ), 
        by = .(pheno_label, probeID)            
    ][,
        median_rank_normalised_reranked := frank(median_rank_normalised, ties.method = "average"), 
        by = .(pheno_label)
    ][
        order(pheno_label, median_rank_normalised_reranked)
    ]

# Plots to check rankings make sense ---------------------------------------------------------------------------------------------------

# make scatter plot comparing median rank to median signed log10 p value
p <- trait_lvl_genes_ranked_by_pSMR_multi %>%
    ggplot(aes(x = median_signed_log10_p_SMR_multi, y = median_rank)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm") +
    facet_wrap(~pheno_label, scales = "free") +
    theme_bw()

p %>% ggsave(filename = "fig_pub/rank_aggregation/median_rank_vs_median_signed_log10_p_SMR_multi_scatterplot.png", width = 20, height = 20)

# normalised rank better reflects the median signed log10 p value
p <- trait_lvl_genes_ranked_by_pSMR_multi %>%
    ggplot(aes(x = median_signed_log10_p_SMR_multi, y = median_rank_normalised)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm") +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~pheno_label, scales = "free") +
    theme_bw()

p %>% ggsave(filename = "fig_pub/rank_aggregation/median_rank_normalised_vs_median_signed_log10_p_SMR_multi_scatterplot.png", width = 20, height = 20)

# -------------------------------------------------------------------------------------------------------------------------------------

phenotypes <- unique(df_msmr_strict$pheno_label)

# get combinations of phenotypes
pheno_combos <- crossing(
    p1 = as.character(phenotypes),
    p2 = as.character(phenotypes)
) %>%
    filter(p1 < p2)

# Calculate spearman rho for all phenotype combinations
spearmancorr_results_summary <- 1:nrow(pheno_combos) %>% 
    map(\(combo_idx) {
        # print(combo_idx)
        p1  <- pheno_combos[combo_idx, 1] %>% pull()
        p2  <- pheno_combos[combo_idx, 2] %>% pull()
        calc_spearman(p1, p2) %>% return()

    }) %>% 
    bind_rows()

spearmancorr_results_summary %>%
    write_tsv(glue("results/rrho/spearman_corr_all_trait_combos_strictmr.tsv"))

# ----------------------------------------------------------------------------------------------------------------------------------
# Range of Spearman correlations for disease traits
# ----------------------------------------------------------------------------------------------------------------------------------

# Calculate stats for the manuscript

trait_meta <- fread("resources/metadata/trait_metadata_n.tsv")
disease_trait_labels <- trait_meta[supercategory == "disease", label]

disease_corrs <- spearmancorr_results_summary %>%
    filter(p1 %in% disease_trait_labels & p2 %in% disease_trait_labels) %>% 
    arrange(spearman_corr)

cat(glue("Disease trait Spearman correlations: min = {min(disease_corrs$spearman_corr)}, max = {max(disease_corrs$spearman_corr)}"), "\n")


# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------




# TODO:
# get hypermat p value BY corrected 
# Is this just the max p value in the BY corrected matrix?

# ----------------------------------------------------------------------------------------------------------------------------------
# Test specific trait pairs
# ----------------------------------------------------------------------------------------------------------------------------------

# p1 <- "Estimated heel BMD"
# p2 <- "Femur neck BMD"

# rrho_out <- calc_rrho_spearman(p1, p2)

# p1 <- "Estimated heel BMD"
# p2 <- "Red blood cell count"

# rrho_obj <- calc_rrho(p1, p2)

# p1 <- "Type 2 diabetes"
# p2 <- "Lung adenocarcinoma"

# rrho_obj <- calc_rrho(p1, p2)

# p1 <- "Bipolar disorder"
# p2 <- "Bipolar I disorder"

# rrho_obj <- calc_rrho(p1, p2)

# p1 <- "ADHD"
# p2 <- "Autism"

# rrho_obj <- calc_rrho(p1, p2)


# Archive
# TEST_USE <- 'two.sided' # options: 'enrichment', 'two.sided'

# calc_rrho <- function(p1, p2, alternative='enrichment') {
#     print(glue("Calculating RRHO for {p1} vs {p2}..."))

#     df1 <- trait_lvl_genes_ranked_by_pSMR_multi %>%
#         filter(pheno_label == p1) %>%
#         select(probeID, median_rank_normalised_reranked)
    
#     df2 <- trait_lvl_genes_ranked_by_pSMR_multi %>%
#         filter(pheno_label == p2) %>%
#         select(probeID, median_rank_normalised_reranked)

#     # Use genes that were tested in both traits
#     intersecting_genes <- intersect(df1$probeID, df2$probeID)
#     df1 <- df1 %>% filter(probeID %in% intersecting_genes)
#     print(head(df1))

#     df2 <- df2 %>% filter(probeID %in% intersecting_genes)
#     print(head(df2))

#     print(glue("RRHO2 on {length(intersecting_genes)} genes that were tested in both traits."))

#     RRHO_obj <- RRHO(
#         df1,
#         df2,
#         BY=TRUE,
#         labels=c(make_clean_names(p1), make_clean_names(p2)),
#         # alternative='enrichment', 
#         alternative=alternative,
#         plots=TRUE,
#         outputdir=glue("/g/data/fy54/analysis/bb3762/repos/tenk10k-causal/fig_pub/rrho_{make_clean_names(alternative)}/"),
#         log10.ind=TRUE)
    
#     max_p_value_raw <- RRHO_obj$hypermat %>% max()
#     max_p_value_by_corrected <- RRHO_obj$hypermat.by %>% max()
#     tibble(
#             p1 = p1,
#             p2 = p2,
#             max_p_value_raw = max_p_value_raw,
#             max_p_value_BY = max_p_value_by_corrected
#         ) %>% return()
# }
# # Calculate RRHO for all phenotype combinations
# rrho_results_summary <- 1:nrow(pheno_combos) %>% 
#     map(\(combo_idx) {
#         print(combo_idx)
#         p1  <- pheno_combos[combo_idx, 1] %>% pull()
#         p2  <- pheno_combos[combo_idx, 2] %>% pull()
#         calc_rrho(p1, p2, alternative=TEST_USE) %>% return()

#     }) %>% 
#     bind_rows()

# rrho_results_summary %>%
#     write_tsv(glue("results/rrho/rrho_pvalues_all_trait_combos_{TEST_USE}.tsv"))    
