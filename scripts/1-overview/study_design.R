# get numbers for manuscript writing

source("scripts/preprocess_strict.R")

library(patchwork)
library(ragg)
library(paletteer)
library(geomtextpath)

# n genes by cell type
df_gene <- df_msmr_tenk10k %>% 
  group_by(cell_type, major_cell_type, probeID, Gene, gene_type) %>% 
  summarise(mr = max(mr), .groups = "drop") %>%
  mutate(gene_cat = ifelse(
    mr == 1, "≥1 trait", "No association"))

df_overall <- df_gene %>% 
  group_by(probeID) %>% 
  slice_max(mr, n = 1, with_ties = FALSE) %>% 
  mutate(cell_type = "Overall", 
         major_cell_type = "Overall")

pals <- paletteer::paletteer_d("MexBrewer::Frida")
(p_gene <- bind_rows(df_gene, df_overall) %>%
  mutate(cell_type = factor(cell_type, c("Overall", levels(df_gene$cell_type))),
         major_cell_type = factor(major_cell_type, c("Overall", levels(df_cell_map$major_cell_type)))) %>% 
  ggplot(aes(x = cell_type, group = fct_rev(gene_cat))) +
  theme_classic() +
  geom_bar(aes(fill = fct_rev(gene_cat))) +
  facet_grid(cols = vars(major_cell_type), scales = "free_x", space = "free_x") +
  # geom_text(aes(label = after_stat(count), group = cell_type,
  #               y = after_stat(count)), vjust = 0,
  #           stat = "count", size = 3) +
  # geom_label(aes(label = after_stat(count)), stat = "count", 
  #            fill = "white", linewidth = 0, alpha = 0.5,
  #            position = position_stack(vjust = 0.5), size = 3) +
  labs(y = "N tested genes", x = NULL) +
  scale_fill_manual(values = pals[c(2,4)], name = "MR association") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(clip = "off") +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 7),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        strip.clip = "off",
        legend.key.size = unit(0.75, "lines"),
        legend.key.spacing.y = unit(0.25, "lines"),
        panel.grid.major.y = element_line(color = "gray90"),
        legend.position = "right")
)

df_traits <- df_msmr_tenk10k[,
  .(n_tested = .SD[,uniqueN(probeID)],
    n_mr = .SD[mr == TRUE, uniqueN(probeID)],
    n_mrsens = .SD[mr_sens == TRUE, uniqueN(probeID)],
    n_mrcoloc = .SD[mr_coloc == TRUE, uniqueN(probeID)],
    n_mrsenscoloc = .SD[mr_sens_coloc == TRUE, uniqueN(probeID)],
    n_magma = .SD[magma_gene == TRUE, uniqueN(probeID)]),
  by = phenotype] |> 
  inner_join(df_trait_map, by = c(phenotype = "trait_id")) |> 
  mutate(pheno_cat = factor(cat_rev, cat_order)) %>% 
  mutate(pheno_label = fct_reorder(label, n_mr), .by = pheno_cat) %>% 
  mutate(scaled_n_eff = scales::rescale(n_eff, to = c(7e3, 1e4)))

(p_traits <- df_traits %>% 
  ggplot(aes(x = pheno_label, fill = pheno_cat)) +
  theme_classic() +
  # geom_point() +
  geom_segment(aes(xend = pheno_label, y = n_mr, yend = -Inf),
               color = "gray90", linetype = "solid") +
  # geom_textline(aes(y = n_mr, label = pheno_cat, group = pheno_cat),
  #               hjust = 0.4, vjust = -0.8, size= 9 /.pt, fontface = "bold") +
  geom_line(aes(y = n_mr, group = pheno_cat)) +
  geom_point(aes(y = n_mr, size = n_eff),
             shape = "circle filled") +
  labs(y = "Number of MR genes", x = NULL) +
  coord_cartesian(clip = "off") +
  # geom_point(aes(y = scaled_n_eff)) +
  facet_grid(cols = vars(pheno_cat),
    labeller = label_wrap_gen(width = 25),
    space = "free", scales = "free") +
  scale_fill_paletteer_d(
    "MetBrewer::Hiroshige", name = NULL,
    # theme = theme(legend.justification = c(0, 0.5)),
    guide = guide_legend(nrow = 2,
      override.aes = list(shape = "square filled", size = 3))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_size_continuous(
    labels = scales::number_format(scale = 1e-3, suffix = "k"),
    name = "Effective GWAS sample size",
    guide = guide_legend(nrow = 1, order = 2, theme = theme(legend.title.position = "left"), override.aes = list(
      fill = "gray", color = "black"))
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 8),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8, margin = margin(l = -1)),
        strip.clip = "off",
        panel.grid.major.y = element_line(color = "gray90"),
        legend.position = "inside",
        legend.key.spacing.y = unit(-0.2, "lines"),
        legend.spacing.y = unit(-0.5, "lines"),
        # legend.box.margin = margin(t = -0.5, b= -0.5),
        plot.margin = margin(),
        legend.position.inside = c(1,1.065),
        legend.justification = c(1,1))
)

infographic <- png::readPNG("figures/study_design/infographics_v3.png", native = TRUE)

plots <- (wrap_elements(full = infographic) /
          wrap_elements(full = p_gene) / 
          wrap_elements(full = p_traits)) +
  plot_layout(heights = c(0.6, 0.2, 0.45)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0, 0.97),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(plots, 
       filename = "figures/study_design/study_design_summary.png", 
       width = 9.5, height = 12, units = "in", dpi = 300, 
       scaling = 1.05,
       device = agg_png)


ggsave(plots, 
       filename = "figures/publication_pdf/Fig1 - Study Design.pdf", 
       width = 9.5/1.05, height = 12/1.05, units = "in", dpi = 320, 
       # scaling = 1.05,
       device = cairo_pdf)

# get numbers
df_traits[supercategory=="disease",.N]
df_gene[,n_distinct(Gene), by = gene_cat]



# write supplementary table
# googlesheets4::gs4_auth()
source("scripts/util/helper.R")

df_celltype <- df_gene %>% 
  group_by(gene_cat, cell_type) %>% 
  tally() %>% 
  ungroup() %>% 
  pivot_wider(names_from = gene_cat, values_from = n, values_fill = 0) %>% 
  mutate(n_gene = `Protein-coding` + `Non-coding`) %>% 
  left_join(df_cell_map, by = "cell_type")

write_gs(df_celltype, "celltype_summary", 1)

df_traits %>% 
  rename(trait_id = phenotype) %>%
  arrange(supercategory, pheno_cat, study, pheno_label) %>%
  mutate(across(starts_with("n_"), as.numeric)) %>% 
  write_gs("trait_summary", 2)


