# MR-coloc overlap: bar chart by phenotype category, and line chart across HEIDI/LFDR thresholds
# combines analyses originally written by Rika and Jayden

library(tidyverse)
library(data.table)
library(ragg)
library(scales)
library(patchwork)

source("scripts/0-preprocess/preprocess_results.R")

# bar chart: average MR-coloc overlap by phenotype category
all_coloc_summary <- df_msmr_tenk10k %>%
  filter(mr == TRUE) %>%
  group_by(phenotype, pheno_cat) %>%
  summarise(
    n_mr = sum(mr, na.rm = TRUE),
    # single OR multi-variant coloc, as defined in the paper
    n_coloc_pph4 = sum(coloc, na.rm = TRUE),
    # n_coloc_pph4 = sum(coloc_pph4 > 0.8, na.rm = TRUE), # single-variant coloc only
    n_coloc_pph3_pph4 = sum(coloc_pph3_h4 > 0.8, na.rm = TRUE),
    # n_mvcoloc_pph4 = sum(mvcoloc_pph4 > 0.8, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # proportion of MR signals overlapping with coloc, for each n_ column except n_mr
  mutate(
    across(
      .cols = c(n_coloc_pph4, n_coloc_pph3_pph4),
      .fns = ~ round(.x / n_mr, digits = 2),
      .names = "{.col}_prop"
    )
  ) %>%
  mutate(sens_inclusive = FALSE)

all_coloc_summary_pheno <- all_coloc_summary %>%
  group_by(pheno_cat) %>%
  summarise(
    average_overlap_pph4 = mean(n_coloc_pph4_prop),
    n = n(),
    n_mean = average_overlap_pph4 * n,
    sd_overlap_pph4 = sd(n_coloc_pph4_prop),
    average_overlap_pph3_pph4 = mean(n_coloc_pph3_pph4_prop),
    sd_overlap_pph3_pph4 = sd(n_coloc_pph3_pph4_prop)
  ) %>%
  mutate(sens_inclusive = FALSE)

# same summary restricted to the sensitivity-analysis MR set
all_coloc_summary_sens <- df_msmr_tenk10k %>%
  filter(mr_sens == TRUE) %>%
  group_by(phenotype, pheno_cat) %>%
  summarise(
    n_mr = sum(mr, na.rm = TRUE),
    n_coloc_pph4 = sum(coloc, na.rm = TRUE),
    n_coloc_pph3_pph4 = sum(coloc_pph3_h4 > 0.8, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      .cols = c(n_coloc_pph4, n_coloc_pph3_pph4),
      .fns = ~ round(.x / n_mr, digits = 2),
      .names = "{.col}_prop"
    )
  ) %>%
  mutate(sens_inclusive = TRUE)

all_coloc_summary_pheno_sens <- all_coloc_summary_sens %>%
  group_by(pheno_cat) %>%
  summarise(
    average_overlap_pph4 = mean(n_coloc_pph4_prop),
    n = n(),
    n_mean = average_overlap_pph4 * n,
    sd_overlap_pph4 = sd(n_coloc_pph4_prop),
    average_overlap_pph3_pph4 = mean(n_coloc_pph3_pph4_prop),
    sd_overlap_pph3_pph4 = sd(n_coloc_pph3_pph4_prop)
  ) %>%
  mutate(sens_inclusive = TRUE)

summary_trait_level <- bind_rows(all_coloc_summary, all_coloc_summary_sens)
summary_phenocat_level <- bind_rows(all_coloc_summary_pheno, all_coloc_summary_pheno_sens)

plot_data_average <- summary_phenocat_level %>%
  mutate(pheno_cat = factor(pheno_cat, levels = unique(pheno_cat))) %>%
  select(pheno_cat, sens_inclusive,
    `H4` = average_overlap_pph4,
    `H3 + H4` = average_overlap_pph3_pph4
  ) %>%
  pivot_longer(
    c(`H4`, `H3 + H4`),
    names_to = "result_type",
    values_to = "proportion"
  ) %>%
  mutate(
    result_type = factor(result_type, levels = c("H4", "H3 + H4")),
    background = factor(
      ifelse(sens_inclusive, "MR + Sensitivity", "MR"),
      levels = c("MR", "MR + Sensitivity")
    )
  )

(p_avg_overlap <- plot_data_average %>%
  filter(result_type == "H4") %>%
  ggplot(aes(x = pheno_cat, y = proportion, fill = background)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_fill_manual(
    values = c("MR" = "#2C7FB8", "MR + Sensitivity" = "#F0A202")
  ) +
  labs(x = NULL, y = "% overlap", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.justification = c(0, 1),
    strip.background = element_rect(fill = "grey90", colour = NA)
  )
)

# ggsave("resources/crohns_case_study/revision/mr_coloc_average_overlap_across_phenocat.png", p_avg_overlap)
# ggsave("resources/crohns_case_study/revision/mr_coloc_average_overlap_across_phenocat_singleonly.png", p_avg_overlap)

# table of average overlap with and without the sensitivity set, by phenotype category
df_coloc_summary_wide <- full_join(all_coloc_summary_pheno, all_coloc_summary_pheno_sens, suffix = c("_mr", "_mr_sens")) %>%
  select(pheno_cat, average_overlap_pph4, sens_inclusive) %>%
  pivot_wider(
    id_cols = pheno_cat,
    names_from = sens_inclusive,
    values_from = c(average_overlap_pph4),
    names_glue = "{ifelse(sens_inclusive, 'mr_sens', 'mr')}"
  )

traits <- all_coloc_summary %>%
  select(phenotype, pheno_cat) %>%
  group_by(pheno_cat) %>%
  summarise(n_traits = n())

df_coloc_summary_wide <- df_coloc_summary_wide %>% left_join(traits, by = "pheno_cat")

# line chart: overlap between MR and coloc support across HEIDI/LFDR thresholds
df_coloc <- df_msmr_tenk10k[, c("cell_type", "phenotype", "probeID", "Gene", "coloc_pph3_h4", "coloc_pph4", "mvcoloc_pph4")]
df_coloc[, coloc_pph3 := coloc_pph3_h4 - coloc_pph4]

df_msmr <- df_msmr_tenk10k[, c("cell_type", "phenotype", "probeID", "Gene", "p_HEIDI", "lfdr_msmr_pheno")]
df_msmr[, neg_log10_q_SMR := -log10(lfdr_msmr_pheno)]

merged_coloc_msmr <- merge(df_coloc, df_msmr, by = c("cell_type", "phenotype", "Gene"))

# parameter ranges
heidi_thresholds <- c(0, 0.01, 0.05, 0.1)
qval_thresholds <- c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5)

# overlap percentages for all HEIDI x LFDR threshold combinations
results_list <- list()

for (heidi in heidi_thresholds) {
  percentages <- sapply(qval_thresholds, function(qval) {
    temp_data <- merged_coloc_msmr %>%
      mutate(
        coloc_sig = ifelse((coloc_pph4 > 0.8) %in% TRUE | (mvcoloc_pph4 > 0.8) %in% TRUE, 1, 0),
        SMR_sig = ifelse(is.na(p_HEIDI), lfdr_msmr_pheno < qval, lfdr_msmr_pheno < qval & p_HEIDI >= heidi)
      )

    n_smr_sig <- sum(temp_data$SMR_sig == 1)
    if (n_smr_sig == 0) return(NA)

    sum(temp_data$SMR_sig == 1 & temp_data$coloc_sig == 1) / n_smr_sig * 100
  })

  percentages_h3 <- sapply(qval_thresholds, function(qval) {
    temp_data <- merged_coloc_msmr %>%
      mutate(
        coloc_h3_sig = ifelse((coloc_pph3 > 0.8) %in% TRUE, 1, 0),
        SMR_sig = ifelse(is.na(p_HEIDI), lfdr_msmr_pheno < qval, lfdr_msmr_pheno < qval & p_HEIDI >= heidi)
      )

    n_smr_sig <- sum(temp_data$SMR_sig == 1)
    if (n_smr_sig == 0) return(NA)

    sum(temp_data$SMR_sig == 1 & temp_data$coloc_h3_sig == 1) / n_smr_sig * 100
  })

  results_list[[as.character(heidi)]] <- data.frame(
    qval = qval_thresholds,
    neg_log_qval = -log10(qval_thresholds),
    percentage = percentages,
    percentage_h3 = percentages_h3,
    heidi_threshold = as.factor(heidi)
  )
}

plot_data <- do.call(rbind, results_list)

(p_h4 <- ggplot(plot_data, aes(x = neg_log_qval, y = percentage, color = heidi_threshold)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("0" = "#984EA3", "0.01" = "#E41A1C", "0.05" = "#377EB8", "0.1" = "#4DAF4A"),
    labels = c(
      bquote(italic(P)[HEIDI] >= 0),
      bquote(italic(P)[HEIDI] >= 0.01),
      bquote(italic(P)[HEIDI] >= 0.05),
      bquote(italic(P)[HEIDI] >= 0.1)
    )
  ) +
  labs(
    x = bquote(-log[10] ~ LFDR[MR] ~ "threshold"),
    y = "% overlap with coloc PP H4 ≥ 0.8",
    color = "HEIDI threshold"
  ) +
  theme_bw() +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  theme(legend.position = "right")
)

(p_h3 <- ggplot(plot_data, aes(x = neg_log_qval, y = percentage_h3, color = heidi_threshold)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("0" = "#984EA3", "0.01" = "#E41A1C", "0.05" = "#377EB8", "0.1" = "#4DAF4A"),
    labels = c(
      bquote(italic(P)[HEIDI] >= 0),
      bquote(italic(P)[HEIDI] >= 0.01),
      bquote(italic(P)[HEIDI] >= 0.05),
      bquote(italic(P)[HEIDI] >= 0.1)
    )
  ) +
  labs(
    x = bquote(-log[10] ~ LFDR[MR] ~ "threshold"),
    y = "% overlap with coloc PP H3 ≥ 0.8",
    color = "HEIDI threshold"
  ) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  theme_bw() +
  theme(legend.position = "right")
)

p_coloc <- ((p_h4 + p_h3) &
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )) +
  plot_layout(guide = "collect")

plots <- (wrap_elements(full = p_avg_overlap) / p_coloc &
  theme(
    panel.grid.minor = element_blank(),
    plot.tag = element_text(face = "bold")
  )) +
  plot_layout(heights = c(5, 4)) +
  plot_annotation(tag_levels = "a")

ggsave("figures/v3/supp/mr_coloc_overlap.png", plots,
  scaling = 1, device = agg_png, height = 8, width = 8
)
