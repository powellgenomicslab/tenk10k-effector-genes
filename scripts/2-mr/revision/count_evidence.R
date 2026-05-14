
# source("scripts/preprocess_strict.R")
suppressPackageStartupMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(arrow)
  library(ragg)
  library(scales)
  library(paletteer)
  library(patchwork)
  library(ggforce)
}))

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v2.parquet.gz")
evidence_criteria <- list(
  mr = expression(lfdr_msmr_pheno < 0.05),
  sensitivity = expression(p_HEIDI >= 0.05 | phet_ivw >= 0.05 | psigmay_mrlink2 >= 0.05),
  coloc = expression(coloc_pph4 >= 0.8 | mvcoloc_pph4 >= 0.8),
  mr_sens = expression(mr & sensitivity),
  mr_coloc = expression(mr & coloc),
  mr_sens_coloc = expression(mr & sensitivity & coloc)
)

# Create columns for evidence criteria and count
for (e in names(evidence_criteria)) {
  df_msmr_tenk10k[, (e) := eval(evidence_criteria[[e]])]
}

df_msmr_tenk10k[, `:=`(
    max_evidence = case_when(
      mr_sens_coloc ~ "mr_sens_coloc",
      mr_coloc      ~ "mr_coloc",
      mr_sens       ~ "mr_sens",
      mr            ~ "mr"
    ) |> factor(levels = c("mr", "mr_sens", "mr_coloc", "mr_sens_coloc"))
  )
]

# add individual evidence expression
evidence <- list(
  mr = expression(lfdr_msmr_pheno < 0.05),
  heidi = expression(p_HEIDI >= 0.05),
  het_ivw = expression(phet_ivw >= 0.05),
  pleio_mrlink2 = expression(psigmay_mrlink2 >= 0.05),
  coloc_single = expression(coloc_pph4 >= 0.8),
  coloc_multi = expression(mvcoloc_pph4 >= 0.8)
)

for (e in names(evidence)) {
  df_msmr_tenk10k[, (e) := eval(evidence[[e]])]
}

evidence_label <- c(
  mr = "MR",
  heidi = "HEIDI test",
  het_ivw = "Cochran's Q test",
  pleio_mrlink2 = "MRLink2 Pleiotropy test",
  coloc_single = "Coloc (single-variant)",
  coloc_multi = "Coloc (multi-variant)"
)

# evidence <- c("mr", "mr_heidi", "coloc_single", "coloc_multi")
max_evidence_label <- c(
  mr = "MR only",
  mr_sens = "MR + Sensitivity",
  mr_coloc = "MR + Coloc",
  mr_sens_coloc = "MR + Sensitivity + Coloc"
)

df_evidence_category <- list(
  MR = c("mr"),
  Sensitivity = c("heidi", "het_ivw", "pleio_mrlink2"),
  Coloc = c("coloc_single", "coloc_multi")
) |> 
  enframe(name = "evidence_category", value = "evidence") |> 
  unnest_longer(evidence) |> setDT()

df_evidence_summary <- df_msmr_tenk10k |> 
  filter(mr) |> 
  select(biosample, phenotype, probeID, all_of(names(evidence)), max_evidence) |> 
  mutate(across(all_of(names(evidence)), ~ifelse(is.na(.x), FALSE, .x) |> as.numeric())) |>
  unite("set", all_of(names(evidence)), sep = "", remove = FALSE) |> 
  setDT()

df_set_tally <- df_evidence_summary |> 
  group_by(set, across(all_of(names(evidence))), max_evidence) |> 
  summarise(N = n(), .groups = "drop") |> 
  arrange(desc(N))

set_order <- unique(df_set_tally$set)

evidence_color <- c(scales::viridis_pal()(4),
"#4E9B8F",
"#E8A025",
"#7B5EA7"
) |> set_names(c(
  "mr", "mr_sens", "mr_coloc", "mr_sens_coloc",
       "MR", "Sensitivity", "Coloc"
))

# tile for set
pals <- paletteer::paletteer_d("MexBrewer::Frida")
# pals <- paletteer::paletteer_d("ggprism::viridis")
# Compute evidence category blocks dynamically from df_evidence_category
# evidence levels (reversed on y-axis): rev(names(evidence)), so position 1 = coloc_multi, ..., 6 = mr
evidence_levels_rev <- rev(names(evidence))
evidence_cat_blocks <- df_evidence_category |>
  mutate(y_pos = match(evidence, evidence_levels_rev)) |>
  summarise(
    ymin = min(y_pos) - 0.5,
    ymax = max(y_pos) + 0.5,
    label_y = mean(y_pos),
    .by = evidence_category
  ) |>
  mutate(evidence_category = factor(evidence_category, levels = c("MR", "Sensitivity", "Coloc")))

cat_fill_cols <- c(MR = "#4E9B8F40", Sensitivity = "#E8A02540", Coloc = "#7B5EA740")

(p_tile <- df_set_tally |> 
  filter(as.numeric(set) > 0) |> 
  pivot_longer(all_of(names(evidence))) |>
  left_join(df_evidence_category, by = c(name = "evidence")) |>
  filter(value == 1) |>
  group_by(set) |>
  mutate(value = na_if(value, 0)) |>
  mutate(set = factor(set, set_order),
         name = factor(name, levels = rev(names(evidence)))) |>
  ggplot(aes(x = set,
             y = name,
             fill = max_evidence)) +
  geom_rect(
    data = evidence_cat_blocks,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf, fill = evidence_category),
    inherit.aes = FALSE,
    alpha = 0.15,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = evidence_color,
    na.translate = FALSE,
    guide = "none"
  ) +
  geom_line(aes(group = set), color = "black", linewidth = 0.5) +
  geom_point(shape = "circle filled", size = 2.5, stroke = 0.5) +
  annotate("label",
           x = -1,
           y = evidence_cat_blocks$label_y,
           label = evidence_cat_blocks$evidence_category,
           hjust = 0.5, angle = 90, vjust = 1,
           fill = "white",
           fontface = "bold",
           color = c("#2E7B70", "#A06A10", "#5A3E87")) +
  theme_void(base_size = 8, base_family = "Helvetica") +
  scale_y_discrete(
    breaks = names(evidence),
    labels = evidence_label,
    expand = expansion(add = 0.5)
  ) +
  scale_x_discrete(expand = expansion(add = c(2.5, 0.5))) +
  theme(panel.grid = element_line(color = "grey90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(hjust = 1, margin = margin(r = 4)),
        plot.margin = margin(l = 40)))

p_set_tally <- df_set_tally |> 
  filter(as.numeric(set) > 0) |> 
  ggplot(aes(x = fct_reorder(set, -N), y = N)) +
  geom_col(aes(fill = max_evidence), color = "black", width = 0.65) +
  geom_text(aes(label = number(N, big.mark = ",", prefix = " ")), nudge_y = 200, angle = 90, hjust = 0, vjust = 0.5) +
  scale_fill_viridis_d(na.translate = FALSE, guide = "none") +
  scale_y_continuous(expand = expansion(c(0, 0.5))) +
  coord_cartesian(clip = "off") +
  layer_scales(p_tile)$x +
  theme_void(base_size = 8, base_family = "Helvetica")

df_evidence_tally <- df_set_tally |> 
  pivot_longer(all_of(names(evidence))) |>
  filter(value == 1) |>
  group_by(name) |> 
  summarise(n_sig = sum(N))

p_evidence_tally <- df_evidence_tally |>
  mutate(name = factor(name, levels = rev(names(evidence)))) |>
  ggplot(aes(y = name, x = n_sig)) +
  geom_col(fill = "gray", color = "black", width = 0.5) +
  scale_x_continuous(expand = expansion(c(0, 1.5))) +
  coord_cartesian(clip = "off") +
  layer_scales(p_tile)$y +
  theme_void(base_size = 8, base_family = "Helvetica") +
  # add label
  geom_text(aes(label = number(n_sig, big.mark = ",")), hjust = -0.1)

# p evidence count tally
(p_evidence_count_tally <- df_evidence_summary |> 
  group_by(max_evidence) |> 
  tally() |> 
  ggplot(aes(y = factor(max_evidence), x = n, fill = factor(max_evidence))) +
  geom_col(color = "black", width = 0.8) +
  scale_fill_viridis_d(guide = "none") +
  scale_x_continuous(expand = expansion(c(0, 0.3))) +
  scale_y_discrete(expand = expansion(c(0, 0)), breaks = names(max_evidence_label), labels = max_evidence_label) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 8, base_family = "Helvetica") +
  # add label
  geom_text(aes(label = number(n, big.mark = ",", prefix = " ")), vjust = 0.5, hjust = 0) +
  labs(y = NULL) +
  theme(axis.text.y = element_text(hjust = 1, margin = margin(r = 4)),
        axis.title.y = element_text(angle = 90, margin = margin(r = 6)))
)
# combine plots
p_combined <- p_set_tally + plot_spacer() + p_tile + p_evidence_tally + plot_layout(ncol = 2, widths = c(3,1)) &
    theme(plot.margin = margin())

p_combined <- wrap_elements(full = p_combined) +
  inset_element(p_evidence_count_tally, left = 0.35, bottom = 0.65, right = 0.85, top = 0.9, align_to = "full")

ggsave("figures/strict/evidence_tally_v2.png",
       p_combined,
       device = ragg::agg_png, scaling = 0.8,
       width = 8, height = 4, units = "in", res = 300)

saveRDS(p_combined, "figures/strict/evidence_tally_plot.rds")

# Individual trait breakdown of max evidence

two_color_pal <- paletteer::paletteer_d("RColorBrewer::RdBu")[c(3,9)] %>% 
  set_names(c("Positive", "Negative"))

df_evidence_by_cat <- df_msmr_tenk10k[mr == TRUE] %>%
  count(cell_type, pheno_cat, supercategory, max_evidence) %>%
  mutate(cat_totals = sum(n), .by = c(pheno_cat, max_evidence)) |> 
  mutate(max_evidence = factor(max_evidence, levels = rev(names(max_evidence_label))))

(p_evidence_by_cat <- df_evidence_by_cat %>%
  ggplot(aes(x = cell_type, y = n, fill = max_evidence)) +
  facet_wrap(~pheno_cat, ncol = 2, scales = "free_y") +
  geom_col() +
  labs(x = NULL, y = "N gene-trait associations") +
  scale_y_continuous(expand = expansion(0),
                     labels = ~ifelse(.x < 1e4, .x, label_number(1, scale = 1e-3, suffix="k")(.x))) +
  scale_fill_manual(values = evidence_color[names(max_evidence_label)],
                    labels = max_evidence_label,
                    breaks = names(max_evidence_label),
                    name = NULL) +
  # scale_alpha_manual(values = setNames(c(0.5, 1), c("Negative", "Positive"))) +
  guides(y = guide_axis(cap = "both"),
         fill = guide_legend(theme = theme(legend.key.size = unit(0.5, "lines")),
                             direction = "horizontal")) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 7),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(colour = "white", fill = "white"),
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "bottom",
    legend.text = element_text(margin = margin(l = 2)),
    # legend.position.inside = c(0.5, 1),
    # legend.justification.inside = c(1, 0),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 9, hjust = 0),
  )
)
ggsave("figures/strict/supp/evidence_by_category.png",
       p_evidence_by_cat,
       device = ragg::agg_png, scaling = 1,
       width = 8, height = 6, units = "in", res = 300)
