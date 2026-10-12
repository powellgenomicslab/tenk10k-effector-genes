
#################################################################################################################
# Interactive job information 
# qsub -I -q normal -P ei56 -l ncpus=2,storage=gdata/fy54+gdata/ei56,mem=32GB -l jobfs=100GB
# cd /g/data/fy54/rt3501/repos/tenk10k-causal/
# conda activate tidyverse
# R
#########################################################################################

library(tidyverse) 
library(ggplot2)
library(ggrepel) 
library(patchwork)

# all_mr <- readRDS("resources/crohns_case_study/postprocess/tenk_alltraits_raw.RDS") %>%
#   filter(pheno_cat == "Immune")
all_mr <- readRDS("resources/crohns_case_study/postprocess/tenk_alltraits_raw.RDS") 
#all_coloc <- readRDS("resources/crohns_case_study/postprocess/df_coloc.RDS") 
#all_mvcoloc <- readRDS("resources/crohns_case_study/postprocess/df_mvcoloc.RDS") 

all_coloc_summary <- all_mr %>% 
  # rename(coloc_pph3 = PP.H3.abf, 
  #        coloc_pph4 = PP.H4.abf) %>% 
  # group_by(phenotype, cell_type) %>% 
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
all_coloc_summary_sens <- all_mr %>% 
  # rename(coloc_pph3 = PP.H3.abf, 
  #        coloc_pph4 = PP.H4.abf) %>% 
  # group_by(phenotype, cell_type) %>% 
  filter(mr_sens == TRUE) %>% ) %>% ################################  NOTE
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

p1 <- ggplot(plot_data_average %>% filter(result_type == "H4"),
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
  #labs(x = "Phenotype Category", y = "Proportion of MR+COLOC Overlap", fill = NULL, title = "Average overlap of MR+COLOC Results By Trait Category (single OR multi-variant)") +
  labs(x = "Phenotype Category", y = "Proportion of MR+COLOC Overlap", fill = NULL, title = "Average overlap of MR+COLOC Results By Trait Category (single variant only)") +
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
df
