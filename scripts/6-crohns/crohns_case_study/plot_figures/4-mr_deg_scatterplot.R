# Purpose: Scatter plot of intersecting MR-DEG genes 

library(tidyverse)
library(ggplot2)
library(ggrepel) 
library(patchwork)

mr_deg <- readRDS("resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS")

plot(mr_deg$b_SMR, mr_deg$`Discrete DE coefficients`)

# cell_map_revised = read_tsv("resources/crohns_case_study/deg/cell_map_revised.tsv")
# colours = cell_map_revised %>% 
#   group_by(revision_major_cell_type) %>% 
#   slice_sample(n = 1) %>% 
#   select(revision_major_cell_type, color) %>% 
#   filter(revision_major_cell_type %in% mr_deg$major_cell_type) %>% 
#   deframe()

darker_colors <- c("#ba8e23", "#674FA3", "#1976D2", "#65350F", "#D84315", "#1B5E20", "#FBC02D", "#4A001F") 


scatterplot <- ggplot(mr_deg, aes(x = p_transform_mr, y = p_transform_deg, colour = major_cell_type, label = Gene)) +
  #scatterplot <- ggplot(mr_deg, aes(x = b_SMR, y = `Discrete DE coefficients`, colour = major_cell_type, label = Gene)) +
  facet_wrap(~major_cell_type, nrow = 4) +
  scale_colour_manual(values = darker_colors) +
  geom_point(alpha = 0.5) + 
  coord_fixed() +
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(n.breaks = 5) +
  geom_text_repel(size = 5, fontface = "italic", max.overlaps = Inf) + 
  theme_minimal() +
  geom_hline(yintercept = 0, linewidth = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.2, alpha = 0.5) +
  # annotate("text", x = -12, y = 12, label = "MR- DEG+", size = 3) +
  # annotate("text",x = 10, y = -12, label = "MR+ DEG-", size = 3) +
  # annotate("text", x = 10, y = 12, label = "MR+ DEG+", size = 3) +
  # annotate("text", x = -12, y = -12, label = "MR- DEG-", size = 3) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  labs(x = bquote(-log[10]~italic(p) %*% "MR direction of effect"), 
       y = bquote(-log[10]~italic(p) %*% "DEG direction of effect")) + 
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold")) +
  ggtitle("Comparison of the Significance and Direction of Effect between TenK10K MR and DEG")

scatterplot

mr_scale <- 18
deg_scale <- 13

scatterplot_all_in_one <- ggplot(mr_deg, aes(x = p_transform_mr, y = p_transform_deg, colour = major_cell_type, label = Gene)) +
  #scatterplot <- ggplot(mr_deg, aes(x = b_SMR, y = `Discrete DE coefficients`, colour = major_cell_type, label = Gene)) +
  # facet_wrap(~major_cell_type) +
  scale_colour_manual(values = darker_colors) +
  geom_point(alpha = 0.5) + 
  coord_fixed() +
  # scale_y_continuous(limits = c(-deg_scale, deg_scale), n.breaks = 10) +
  # scale_x_continuous(limits = c(-mr_scale, mr_scale), n.breaks = 10) +
  scale_y_continuous(limits = c(-9, 14), n.breaks = 10) +
  scale_x_continuous(limits = c(-17.5, 13), n.breaks = 10) +
  geom_text_repel(key_glyph = "point", size = 3.5, fontface = "italic", max.overlaps = Inf) + 
  theme_minimal() +
  geom_hline(yintercept = 0, linewidth = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.2, alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, 
           alpha = 0.1, fill = "lightgreen") + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "#f9bfbf") + 
  annotate("text", x = -16, y = 14, label = "MR- DEG+", size = 4) +
  annotate("text",x = 11.5, y = -9, label = "MR+ DEG-", size = 4) +
  annotate("text", x = 11.5, y = 14, label = "MR+ DEG+", size = 4) +
  annotate("text", x = -16, y = -9, label = "MR- DEG-", size = 4) +
  labs(x = bquote(-log[10] ~ italic(P)[MR] %*% "direction of effect"), 
       y = bquote(-log[10] ~ italic(P)[DEG] %*% "direction of effect"),
       colour = "Major Cell Type") + 
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) 
  # ggtitle("Comparison of the Significance and Direction of Effect between TenK10K MR and DEG")

scatterplot_all_in_one

save(scatterplot, scatterplot_all_in_one, file = "resources/crohns_case_study/figures/mr_deg_scatterplot_objects.RData")

ggsave("resources/crohns_case_study/figures/mr_deg_facet_annot.png", scatterplot, device = ragg::agg_png(),
       width = 5.0, height = 8, bg = "white", scaling = 0.5, dpi = 300)

ggsave("resources/crohns_case_study/figures/mr_deg.png", scatterplot_all_in_one, device = ragg::agg_png(),
       width = 5.0, height = 4.5, bg = "white", scaling = 0.7, dpi = 300)

