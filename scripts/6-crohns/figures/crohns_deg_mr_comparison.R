# Compare direction of effect between DEG in crohn's and MR
library(data.table)
library(tidyverse)
library(readxl)
library(ragg)
library(paletteer)
library(ggrepel)

# read in pre-processed Crohn's sig results  
df_msmr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% ungroup() %>% as.data.table()
# pre-processed deg results 
df_deg <- readRDS("resources/crohns_case_study/deg/crohns_deg_pre-processed.RDS") %>% ungroup() %>% as.data.table()

# get data for plotting (intersection MR - DEG results)
df_mr_deg <- df_msmr %>%
  filter(phenotype == "crohns") %>%
  select(Gene, cell_type, b_SMR, se_SMR, p_SMR, p_SMR_multi, p_HEIDI, major_cell_type) %>%
  inner_join(df_deg, by = c("Gene", "major_cell_type")) %>%
  mutate(concordant_mr_deg = sign(b_SMR) == sign(`Discrete DE coefficients`)) %>%
  mutate(concordance_group = case_when(
            all(concordant_mr_deg) ~ "All concordant",
            all(!concordant_mr_deg) ~ "All discordant",
            any(concordant_mr_deg) & any(!concordant_mr_deg) ~ "Mixed concordance"
          ), .by = c("Gene"))


# cross plot
df_plot <- bind_rows(
  df_mr_deg[concordance_group != "Mixed concordance", .(top_b_smr = .SD[which.min(p_SMR_multi), b_SMR],
              top_p_smr = min(p_SMR_multi),
              top_b_deg = .SD[which.min(`Discrete DE coefficients p value`),
                              `Discrete DE coefficients`],
              top_p_deg = min(`Discrete DE coefficients p value`)),
          by = c("Gene", "major_cell_type", "concordance_group")],
  df_mr_deg[concordance_group == "Mixed concordance", .(top_b_smr = .SD[which.min(p_SMR_multi), b_SMR],
                                                        top_p_smr = min(p_SMR_multi),
                                                        top_b_deg = .SD[which.min(`Discrete DE coefficients p value`),
                                                                        `Discrete DE coefficients`],
                                                        top_p_deg = min(`Discrete DE coefficients p value`)),
            by = c("Gene", "major_cell_type", "concordance_group", "concordant_mr_deg")]
)

max_b = 5

gene_highlight <- "TNFRSF18"
(p_cross <- df_mr_deg %>% 
    mutate(b_smr = ifelse(abs(b_SMR) > max_b, sign(b_SMR) * max_b, b_SMR),
           b_deg = ifelse(abs(`Discrete DE coefficients`) > max_b, sign(`Discrete DE coefficients`) * max_b, `Discrete DE coefficients`),
          # gene_label = ifelse(Gene %in% gene_highlight, 
          #                      paste(Gene, cell_type, scRNAseq_cellid, sep = " | "), "")) %>%
          gene_label = pmap(list(Gene, cell_type, scRNAseq_cellid),
                            ~if (..1 %in% gene_highlight) {
                            bquote(italic(.(..1))[.(as.character(..2)) ~ "|" ~ .(..3)])
                           } else "")) %>% 
    
            # , 
            #                    paste(Gene, cell_type, scRNAseq_cellid, sep = " | "), "")) %>%
  # mutate(top_b_smr = ifelse(abs(top_b_smr) > max_b, sign(top_b_smr) * max_b, top_b_smr),
  #        top_b_deg = ifelse(abs(top_b_deg) > max_b, sign(top_b_deg) * max_b, top_b_deg),
  #.       gene_label = ifelse(Gene %in% gene_highlight, Gene, "")) %>%
  # ggplot(aes(y = top_b_smr, x = top_b_deg, color = major_cell_type)) +
  ggplot(aes(y = b_smr, x = b_deg, color = major_cell_type)) +
  geom_point(aes(fill = after_scale(alpha(color,  0.6))),
             size = 2, shape = "circle filled") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~concordance_group, scales = "fixed", nrow = 1,
             axes = "all") +
  labs(x = bquote(beta[DEG]),
       y = bquote(beta[MR])) +
  geom_text_repel(aes(label = gene_label), segment.size = 0.2,
                  data = ~.x[Gene %in% gene_highlight],
                  parse = TRUE,
                   color = "black", xlim = c(0, 5),
                   force = 4, hjust = 0.5,
                   ylim = c(1, 5), seed =99,  size = 8/.pt,
                   max.overlaps = Inf) +
  geom_label(
    aes(label = paste("N genes:", n)),
    linewidth = 0, size = 9/.pt,
    inherit.aes = FALSE,
    data = ~.x[,.(n = n_distinct(Gene)), by = concordance_group],
    x = -max_b, y = max_b,
    hjust = 0, vjust = 1
  ) +
  scale_x_continuous(labels = ~case_when(.x == -max_b ~ paste0("≤", .x),
                                            .x == max_b ~ paste0("≥", .x),
                                            TRUE ~as.character(.x))) +
  scale_y_continuous(labels = ~case_when(.x == -max_b ~ paste0("≤", .x),
                                            .x == max_b ~ paste0("≥", .x),
                                            TRUE ~as.character(.x))) +
  scale_color_paletteer_d("ggthemes::Classic_10", name = "Major cell type") +
  theme_bw() +
  coord_fixed(xlim = c(-1, 1) * max_b,
              ylim = c(-1, 1) * max_b,
              clip = "off") +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
)
ggsave("figures/crohns/crohns_deg_mr_cross_plot.png", p_cross,
       width = 8, height = 4, device = agg_png, scaling = 0.8)



# #add the "tenk10k harmonised" cell types - basically their cell type groups 
# cell_map <- read.delim("resources/metadata/cell_map.tsv")
# major_cell_types_tenk <-  unique(cell_map$major_cell_type)
# 
# deg$major_cell_type <- deg_cell_map$major_cell_type[match(deg$scRNAseq_cellid, deg_cell_map$scRNAseq_cellid)]
# 
# deg$major_cell_type <- as.factor(deg$major_cell_type)
# deg$scRNAseq_cellid <- as.factor(deg$scRNAseq_cellid)
# 
# 
# # create a text file of the common cell types maybe useful later.
# # intersect major cell types from tenk and deg
# common_cell_types <- intersect(cell_map$major_cell_type, deg_cell_map$major_cell_type)
# 
# deg_crohns_supp <- readxl::read_excel("resources/misc/Crohns_Summary.xlsx",
#                                   sheet = "crohns_annotated_results_deg") %>% 
#   select(probeID, cell_type,
#          major_cell_type, Gene, Location, Contrast,
#          scRNAseq_cellid = crohns_dataset_scRNAseq_cell_id,
#          concordant_supp = concordant_DEG_and_MR_direction) %>%
#   # rename(scRNAseq_cellid = crohns_dataset_scRNAseq_cell_id,
#   #        concordant_supp = concordant_DEG_and_MR_direction) %>%
#   filter(n_distinct(concordant_supp) > 1, .by = c(probeID, cell_type,
#            major_cell_type, Gene, Location, Contrast,
#            scRNAseq_cellid)) %>% 
#   arrange(probeID, cell_type,
#           major_cell_type, Gene, Location, Contrast,
#           scRNAseq_cellid)
# 
# deg_crohns <- readRDS("resources/crohns_case_study/deg/MR_and_DEG_matched_major_cell_type_combined_results_colon_and_ti_df_concordance_annotated.RDS") %>% 
#   select(probeID, cell_type,
#          major_cell_type, Gene, Location, Contrast,
#          scRNAseq_cellid,
#          concordant)
# 
# inner_join(deg_crohns_supp, deg_crohns) %>% 
#   filter(concordant_supp != concordant) %>%
#   View()
# 
# 
# 
# 
# anti_join(deg_crohns_supp, deg_crohns,
#           by = c("probeID", "cell_type", "concordant_DEG_and_MR_direction" = "concordant",
#                  "major_cell_type", "Gene", "Location", "Contrast",
#                  "crohns_dataset_scRNAseq_cell_id" = "scRNAseq_cellid"))
# 
# df_plot_deg <- deg_crohns %>% 
#   group_by(Gene, major_cell_type)
# 
# setDT(deg_crohns)
# df_crohns_major <- df_msmr[phenotype == "crohns", .SD[which.min(p_SMR_multi)],
#         by = .(major_cell_type, probeID)]
# 
# deg_crohns[df_crohns_major, `:=`(b_smr = i.b_SMR,
#                   p_msmr = i.p_SMR_multi),
#            on = c("major_cell_type", "probeID")]
# 
# deg_crohns_discrete <- deg_crohns[Discrete.FDR < 0.05] %>% 
#   group_by(Gene, major_cell_type) %>% 
#   slice_min(Discrete.DE.coefficients.p.value, with_ties = FALSE) %>% 
#   mutate(deg_model = "discrete",
#          beta_deg = Discrete.DE.coefficients,
#          p_deg = Discrete.DE.coefficients.p.value)
# deg_crohns_cont <- deg_crohns[Continuous.FDR < 0.05] %>% 
#   group_by(Gene, major_cell_type) %>% 
#   slice_min(Continuous.DE.coefficients.p.value, with_ties = FALSE) %>% 
#   mutate(deg_model = "cont",
#          beta_deg = Continuous.DE.coefficients,
#          p_deg = Continuous.DE.coefficients.p.value)
# 
# deg_crohns_plot <- bind_rows(deg_crohns_discrete, deg_crohns_cont) %>% 
#   ungroup() %>% 
#   select(-cell_type) %>% 
#   distinct() %>% 
#   mutate(prop_concordant = sum(sign(b_smr) ==  sign(beta_deg)) / n(),
#          .by = deg_model)
# 
# deg_crohns_plot %>% 
#   summarise(prop_concordant = sum(sign(b_smr) ==  sign(beta_deg)) / n(),
#             n = n(),
#          .by = c(deg_model, major_cell_type)) %>% 
#   summarise(mean_prop = mean(prop_concordant), .by = deg_model)
#             n = sum(n)
#   pull(prop_concordant) %>%
#   mean()
# deg_crohns_plot %>% 
#   # filter(major_cell_type == "B") %>% 
#   ggplot(aes(y = b_smr, x = beta_deg, color = major_cell_type)) +
#   geom_point() +
#   facet_grid(cols = vars(deg_model),
#              rows = vars(major_cell_type),
#              scales = "free") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   theme_bw() +
#   theme(aspect.ratio = 1)
# 
