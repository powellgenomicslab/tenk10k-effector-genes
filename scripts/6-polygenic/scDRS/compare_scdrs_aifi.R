# compare scDRS between TenK10K and Immune Health Atlas (AIFI) datasets

library(data.table)
library(tidyverse)
library(scales)
library(patchwork)
library(ragg)
library(ggforce)
library(ggthemes)
library(ggridges)
library(ggrepel)
library(readxl)
library(glue)
library(geomtextpath)
library(paletteer)

df_cell_map <- fread("metadata/cell.tsv") %>% 
  mutate(cell_type = factor(cell_type, unique(cell_type)))

df_trait_map <- fread("metadata/trait.tsv") %>% 
  filter(include) %>% 
  setDT()
df_trait_cat <- fread("metadata/trait_category.tsv") %>% 
  setDT() 

df_stats <- list(
  `tenk10k` = "results/aggregate/tenk10k_phase1.scdrs.cell_type_stats.tsv",
  `aifi` = "results/aggregate/immune_health_atlas_annotated.scdrs.cell_type_stats.tsv"
) %>% 
  map_df(fread, .id = "study")

# compare number of cells
p_cell_count <- df_stats |> 
  distinct(study, cell_type, n_cell) |> 
  mutate(study = case_when(
    study == "tenk10k" ~ "TenK10K",
    study == "aifi" ~ "AIFI Immune Health Atlas",
    TRUE ~ study
  ),
         cell_type = str_replace_all(cell_type, "_", " ")) |>
  ggplot(aes(x = fct_reorder(cell_type, n_cell, diff), fill = study, y = n_cell)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  theme_minimal() +
  scale_fill_tableau() +
  labs(x = NULL, y = "Number of sampled cells", fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "top",
        # legend.background = element_rect(fill = "white", color = NA),
        # legend.position.inside = c(0,1),
        legend.justification = c(0,0),
        legend.key.size = unit(0.5, "lines"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5)
      )

# Supplementary Figure 23: cells sampled from AIFI and TenK10K by cell type
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig-s23.png", p_cell_count,
       device = agg_png, scaling = 1,
        width = 8, height = 5, units = "in", res = 300)
ggsave("figures/fig-s23.pdf", p_cell_count,
       device = cairo_pdf,
        width = 8, height = 5, units = "in")


# compare results
df_stats[df_trait_map, `:=`(pheno_label = i.label,
                            pheno_cat = i.cat_rev,
                            pheno_slim = i.include,
                            supercategory = i.supercategory),
         on = c(phenotype = "trait_id")]
df_stats[, pheno_cat := factor(pheno_cat, df_trait_cat$cat_order)]

df_stats_wide <- df_stats %>% 
  filter(complete.cases(.)) %>%
  select(phenotype, cell_type, assoc_mcz, assoc_mcp, n_cell, study) %>%
  pivot_wider(names_from = study, values_from = c(assoc_mcz, assoc_mcp, n_cell), names_sep = ".") %>%
  left_join(df_trait_map %>% select(phenotype = trait_id, phenotype_label = label, cat_rev, supercategory),
            by = "phenotype") 
  # %>% 
  # filter to cell-type significant in TenK10K
  # filter(assoc_mcp.tenk10k < 0.05)
 
# find how many cell types have higher mcz
df_stats |> 
  ggplot(aes(x = study, group = study, y = abs(assoc_mcz))) +
  facet_wrap(~cell_type) +
  geom_boxplot() +
  theme_minimal()
# calculate correlation by phenotype
cor_mcz_pheno <- df_stats_wide %>%
  group_by(phenotype) %>%
  summarise(cor_z = cor(assoc_mcz.tenk10k, assoc_mcz.aifi, method = "pearson", use = "complete.obs"),
            n = n())

cor_mcz_celltype <- df_stats_wide %>%
  group_by(cell_type) %>%
  summarise(cor_z = cor(assoc_mcz.tenk10k, assoc_mcz.aifi, method = "pearson", use = "complete.obs"))

xy_limits <- range(c(df_stats_wide$assoc_mcz.tenk10k, df_stats_wide$assoc_mcz.aifi), na.rm = TRUE)

(
p_correlation <- df_stats_wide %>% 
  ggplot(aes(x = assoc_mcz.tenk10k, y = assoc_mcz.aifi)) +
  geom_point(aes(color = supercategory), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  # geom_smooth(method = "lm", color = "steelblue", linewidth = 0.8, se = FALSE, fullrange = TRUE) +
  facet_wrap(~cell_type, nrow = 5, labeller = labeller(cell_type = ~str_replace_all(.x, "_", " "))) +
  # annotate with correlation
  geom_label(data = cor_mcz_celltype,
             aes(x = xy_limits[1] ,
                 y = xy_limits[2] ,
                 label = paste0("r = ", round(cor_z, 2))),
             linewidth = 0, fill = alpha("white", 0.5),
             label.padding = unit(0.1, "lines"),
             hjust = 0, vjust = 1,
             inherit.aes = FALSE,
             size = 8/.pt,
             color = "black") +
  labs(x = "TenK10K scDRS Z-score",
       y = "AIFI scDRS Z-score",
       title = "Comparison of scDRS MCZ association between TenK10K and AIFI datasets") +
  scale_color_paletteer_d("ggthemes::Tableau_10", name = NULL, labels = str_to_title) +
  theme_bw() +
  coord_fixed(xlim = xy_limits, ylim = xy_limits) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(1,0),
        legend.justification = c(1,0))
)
# save plot
# Supplementary Figure 24: scDRS Z-score, TenK10K vs AIFI
ggsave("figures/fig-s24.png", p_correlation,
       device = agg_png, scaling = 1.25,
        width = 10, height = 10, units = "in", res = 300)
ggsave("figures/fig-s24.pdf", p_correlation,
       device = cairo_pdf, scale = 1/1.25,
        width = 10, height = 10, units = "in")

# correlation by trait
ggplot(cor_mcz_pheno, aes(y = cor_z, x = phenotype)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(y = "Correlation of scDRS MCZ by phenotype",
       x = "Phenotype") +
  theme_bw()

# compare number of cells
df_stats_wide