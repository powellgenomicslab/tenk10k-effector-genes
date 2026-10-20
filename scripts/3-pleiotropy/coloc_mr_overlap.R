
#################################################################################################################
# Mix of Rika's and Jayden script
# #########################################################################################

library(tidyverse) 
library(ggplot2)
library(ggrepel) 
library(patchwork)

source("scripts/0-preprocess/preprocess_results.R")

# BAR PLOT

all_coloc_summary <- df_msmr_tenk10k %>% 
  filter(mr == TRUE) %>% 
  group_by(phenotype, pheno_cat) %>% 
  summarise(
    n_mr = sum(mr, na.rm = TRUE),
    n_coloc_pph4 = sum(coloc, na.rm = TRUE), ################################ NOTE this is single | multi i.e. definition in paper 
    #n_coloc_pph4 =sum(coloc_pph4 >0.8, na.rm = TRUE), ################################ NOTE this is single only 
    n_coloc_pph3_pph4 = sum(coloc_pph3_h4 >0.8, na.rm = TRUE),
    #n_mvcoloc_pph4 = sum(mvcoloc_pph4 > 0.8, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  # Dynamically calculate proportions for all 'n_' columns except n_mr
  mutate(
    across(
      #.cols  = c(n_coloc_pph4, n_coloc_pph3_pph4, n_mvcoloc_pph4), 
      .cols  = c(n_coloc_pph4, n_coloc_pph3_pph4), 
      .fns   = ~ round(.x / n_mr, digits = 2), 
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

################################################ Summary with Sensitivity Analysis 
all_coloc_summary_sens <- df_msmr_tenk10k %>% 
  # rename(coloc_pph3 = PP.H3.abf, 
  #        coloc_pph4 = PP.H4.abf) %>% 
  # group_by(phenotype, cell_type) %>% 
  filter(mr_sens == TRUE) %>% 
  group_by(phenotype, pheno_cat) %>% 
  summarise(
    n_mr = sum(mr, na.rm = TRUE),
    n_coloc_pph4 = sum(coloc, na.rm = TRUE), ################################ NOTE this is single | multi i.e. definition in paper  
    #n_coloc_pph4 =sum(coloc_pph4 >0.8, na.rm = TRUE), ################################ NOTE this is single only 
    n_coloc_pph3_pph4 = sum(coloc_pph3_h4 >0.8, na.rm = TRUE),
    #n_mvcoloc_pph4 = sum(mvcoloc_pph4 > 0.8, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  # Dynamically calculate proportions for all 'n_' columns except n_mr
  mutate(
    across(
      #.cols  = c(n_coloc_pph4, n_coloc_pph3_pph4, n_mvcoloc_pph4), 
      .cols  = c(n_coloc_pph4, n_coloc_pph3_pph4), 
      .fns   = ~ round(.x / n_mr, digits = 2), 
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
################################################ Summary with Sensitivity Analysis 

summary_trait_level <- bind_rows(all_coloc_summary, all_coloc_summary_sens)
summary_phenocat_level <- bind_rows(all_coloc_summary_pheno, all_coloc_summary_pheno_sens)

# calculate delta
df_summary_phenocat_wide <- summary_phenocat_level |> 
  mutate(sens_inclusive = ifelse(sens_inclusive, "MR + Sensitivity", "MR")) |> 
  select(pheno_cat, name = sens_inclusive, value = average_overlap_pph4) |> 
  pivot_wider() |> 
  mutate(delta_average = `MR + Sensitivity` - MR)

plot_data_average <- summary_phenocat_level %>%
  # Order phenotypes by category so each category forms a contiguous block on the axis
  #arrange(pheno_cat, phenotype) %>%
  mutate(pheno_cat = factor(pheno_cat, levels = unique(pheno_cat))) %>%
  select(pheno_cat, sens_inclusive,
         `H4` = average_overlap_pph4,
         `H3 + H4` = average_overlap_pph3_pph4) %>%
  pivot_longer(
    c(`H4`, `H3 + H4`),
    names_to  = "result_type",
    values_to = "proportion"
  ) %>%
  mutate(
    result_type = factor(result_type, levels = c("H4", "H3 + H4")),
    background  = factor(
      ifelse(sens_inclusive, "MR + Sensitivity", "MR"),
      levels = c("MR", "MR + Sensitivity")
    ),
    # Nested axis label: phenotype on top line, category bracketed below
    #x_label = interaction(phenotype, pheno_cat, sep = " | ", lex.order = FALSE)
  )
(
p1 <- ggplot(plot_data_average %>% filter(result_type == "H4"),
       aes(x = pheno_cat, y = proportion, fill = background)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  #facet_wrap(~ result_type, nrow = 1) +
  # scale_x_discrete(guide = guide_axis_nested(delim = "|")) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(seq(0, 1, 0.2))
  ) +
  scale_fill_manual(
    # values = c("MR" = "#4E9B8F", "MR + Sensitivity" = "#E8A025")
    values = c("MR" = "#4E79A7FF", "MR + Sensitivity" = "#F28E2BFF")
  ) +
  #labs(x = "Phenotype Category", y = "Proportion of MR+COLOC Overlap", fill = NULL, title = "Average overlap of MR+COLOC Results By Trait Category (single OR multi-variant)") +
  labs(x = NULL, y = "% overlap", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x           = element_text(angle = 30, hjust = 1, size = 9),
    #ggh4x.axis.nesttext.x = element_text(face = "bold", size = 9),
    #ggh4x.axis.nestline.x = element_line(linewidth = 0.5),
    panel.grid.major.x    = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.spacing.x       = unit(0.5, "lines"),
    legend.position       = "inside",
    legend.position.inside = c(0.01,0.99),
    legend.justification = c(0,1),
    strip.background      = element_rect(fill = "grey90", colour = NA)
  )
)



#ggsave("resources/crohns_case_study/revision/mr_coloc_average_overlap_across_phenocat.png", p1)
#ggsave("resources/crohns_case_study/revision/mr_coloc_average_overlap_across_phenocat_singleonly.png", p1)

########################################################################################################################################

df <- full_join(all_coloc_summary_pheno, all_coloc_summary_pheno_sens, suffix = c("_mr", "_mr_sens")) %>% select(pheno_cat, average_overlap_pph4, sens_inclusive)
df <- df%>%
  pivot_wider(
    id_cols     = pheno_cat,
    names_from  = sens_inclusive,
    values_from = c(average_overlap_pph4),
    names_glue  = "{ifelse(sens_inclusive, 'mr_sens', 'mr')}"
  )
traits <- all_coloc_summary %>% select(phenotype, pheno_cat) %>% group_by(pheno_cat) %>% summarise(n_traits = n())
df <- df %>% left_join(traits, by = "pheno_cat")

# LINE CHART From Jayden

df_coloc <- df_msmr_tenk10k[, c("cell_type", "phenotype", "probeID", "Gene", "coloc_pph3_h4", "coloc_pph4", "mvcoloc_pph4")]
df_coloc$coloc_pph3 <- df_coloc$coloc_pph3_h4 - df_coloc$coloc_pph4

df_msmr <- df_msmr_tenk10k[, c("cell_type", "phenotype", "probeID", "Gene", "p_HEIDI", "lfdr_msmr_pheno")]
df_msmr[, neg_log10_q_SMR := -log10(lfdr_msmr_pheno)]

merged_coloc_msmr <- merge(
  df_coloc,
  df_msmr,
  by.x = c("cell_type", "phenotype", "Gene"),
  by.y = c("cell_type", "phenotype", "Gene")
)

# Define parameter ranges
heidi_thresholds <- c(0, 0.01, 0.05, 0.1)
qval_thresholds <- c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5)

# sequential palette: darker shade = more stringent HEIDI threshold
# (magma runs dark -> light by default, so reverse it: light = 0, dark = 0.1)
heidi_threshold_color <- setNames(
  rev(viridisLite::magma(length(heidi_thresholds), begin = 0.15, end = 0.85)),
  as.character(heidi_thresholds)
)

# Calculate overlap percentages for all combinations
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

# Combine results
plot_data <- do.call(rbind, results_list)

# Create line chart
p_h4 <- ggplot(plot_data, aes(x = neg_log_qval, y = percentage, color = heidi_threshold)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = heidi_threshold_color,
    labels = c(
      bquote(italic(P)[HEIDI] >= 0),
      bquote(italic(P)[HEIDI] >= 0.01),
      bquote(italic(P)[HEIDI] >= 0.05),
      bquote(italic(P)[HEIDI] >= 0.1)
  )) +
  labs(
    x = bquote(-log[10]~LFDR[MR]~"threshold"),
    # title = "Proportion of MR signals\noverlapping with coloc PP H4 ≥ 0.8",
    y = "% overlap with coloc PP H4 ≥ 0.8",
    color = "HEIDI threshold"
  ) +
  theme_bw() +
  scale_y_continuous(labels = ~paste0(.x,"%"))+
  theme(legend.position = "right")

p_h3 <- ggplot(plot_data, aes(x = neg_log_qval, y = percentage_h3, color = heidi_threshold)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = heidi_threshold_color,
    labels = c(
      bquote(italic(P)[HEIDI] >= 0),
      bquote(italic(P)[HEIDI] >= 0.01),
      bquote(italic(P)[HEIDI] >= 0.05),
      bquote(italic(P)[HEIDI] >= 0.1)
  )) +
  labs(
    x = bquote(-log[10]~LFDR[MR]~"threshold"),
    # title = "Proportion of MR signals\noverlapping with coloc PP H3 ≥ 0.8",
    y = "% overlap with coloc PP H3 ≥ 0.8",
    color = "HEIDI threshold"
  ) +
  scale_y_continuous(labels = ~paste0(.x,"%"))+
  theme_bw() +
  theme(legend.position = "right")


p_coloc <- ((p_h4 + p_h3) &
  theme(panel.grid_minor = element_blank(),
        legend.position = "bottom")) +
        # legend.position.inside = c(0.01,0.99),
        # legend.justification = c(0,1)) +
  plot_layout(guide = "collect") 

plots <-   (wrap_elements(full = p1) / p_coloc &
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold"))) +
  plot_layout(heights = c(5,4)) +
  plot_annotation(tag_levels = "a")

# Supplementary Fig. 8a-c
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig-s8.png", plots, scaling = 1,
       device = ragg::agg_png, height = 8, width = 8)
ggsave("figures/fig-s8.pdf", plots,
       device = cairo_pdf, height = 8, width = 8)
