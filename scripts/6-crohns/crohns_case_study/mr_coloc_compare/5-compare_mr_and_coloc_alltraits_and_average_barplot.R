library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)

plot_data <- summary_trait_level %>%
  # Order phenotypes by category so each category forms a contiguous block on the axis
  arrange(pheno_cat, phenotype) %>%
  mutate(phenotype = factor(phenotype, levels = unique(phenotype))) %>%
  select(phenotype, pheno_cat, sens_inclusive,
         `H4` = n_coloc_pph4_prop,
         `H3 + H4` = n_coloc_pph3_pph4_prop) %>%
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
    x_label = interaction(phenotype, pheno_cat, sep = " | ", lex.order = FALSE)
  )

ggplot(plot_data,
       aes(x = x_label, y = proportion, fill = background)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ result_type, nrow = 1) +
  scale_x_discrete(guide = guide_axis_nested(delim = "|")) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(0, 0.5, 1)
  ) +
  scale_fill_manual(
    values = c("MR" = "#2C7FB8", "MR + Sensitivity" = "#F0A202")
  ) +
  labs(x = NULL, y = "Proportion of results", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x           = element_text(angle = 45, hjust = 1, size = 6),
    ggh4x.axis.nesttext.x = element_text(face = "bold", size = 9),
    ggh4x.axis.nestline.x = element_line(linewidth = 0.5),
    panel.grid.major.x    = element_blank(),
    panel.spacing.x       = unit(0.5, "lines"),
    legend.position       = "bottom",
    strip.background      = element_rect(fill = "grey90", colour = NA)
  )

###############################################################################################################################

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

ggplot(plot_data_average %>% filter(result_type == "H4"),
       aes(x = pheno_cat, y = proportion, fill = background)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  #facet_wrap(~ result_type, nrow = 1) +
  scale_x_discrete(guide = guide_axis_nested(delim = "|")) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(seq(0, 1, 0.1))
  ) +
  scale_fill_manual(
    values = c("MR" = "#2C7FB8", "MR + Sensitivity" = "#F0A202")
  ) +
  labs(x = "Phenotype Category", y = "Proportion of MR+COLOC Overlap", fill = NULL, title = "Average overlap of MR+COLOC Results By Trait Category (single-variant only)") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x           = element_text(angle = 45, hjust = 1, size = 10),
    #ggh4x.axis.nesttext.x = element_text(face = "bold", size = 9),
    #ggh4x.axis.nestline.x = element_line(linewidth = 0.5),
    #panel.grid.major.x    = element_blank(),
    #panel.spacing.x       = unit(0.5, "lines"),
    legend.position       = "bottom",
    strip.background      = element_rect(fill = "grey90", colour = NA)
  )

