# Author: Blake Bowen

library(glue)
library(tidyverse)
library(arrow)
library(data.table)

# fixed sort order so trait pairs (p1 < p2) are written the same way regardless of IDE/terminal locale
Sys.setlocale("LC_COLLATE", "C")

dir.create("results/rrho", recursive = TRUE, showWarnings = FALSE)

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
    df2 <- df2 %>% filter(probeID %in% intersecting_genes)

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

df_msmr_strict <- read_parquet("results/preprocessed/tenk10k_phase1.v5.parquet.gz") %>%
    setDT()

# start with the full mr result list - not filtered for sig only
#  use strict msmr results, which are stringently filtered to only include more significant eGenes
pheno_order <- df_msmr_strict[, .N, by = pheno_label][order(-N), pheno_label]
df_msmr_strict[, pheno_label := factor(pheno_label, pheno_order)]
df_msmr_strict <- df_msmr_strict[p_SMR_multi < 0.05,] # filter to keep only nominally significant genes to reduce noise

# ----------------------------------------------------------------------------------------------------------------------------------
# get significant gene sets for each trait
# ----------------------------------------------------------------------------------------------------------------------------------

# all cell types combined
# method for combining gene lists across cell types, to get a single ranked gene list per trait:
# 1) rank genes by signed -log10(p_SMR_multi) for each trait and each cell type, take the median rank across cell types

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

trait_meta <- fread("metadata/trait.tsv")
disease_trait_labels <- trait_meta[supercategory == "disease", label]

disease_corrs <- spearmancorr_results_summary %>%
    filter(p1 %in% disease_trait_labels & p2 %in% disease_trait_labels) %>%
    arrange(spearman_corr)

cat(glue("Disease trait Spearman correlations: min = {min(disease_corrs$spearman_corr)}, max = {max(disease_corrs$spearman_corr)}"), "\n")
