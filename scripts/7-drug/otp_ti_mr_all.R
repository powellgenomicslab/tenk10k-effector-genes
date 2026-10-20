# Extended Data: Open Targets association scores of target-indication (T-I) pairs
# with MR support, by disease. The top 15 targets per disease are labelled and
# points are coloured by maximum clinical development phase.
# Moved from otp_combined.R; the data preparation below is copied from there
# (keep in sync).
source("scripts/0-preprocess/preprocess_results.R")

library(ggrepel)
library(ragg)
library(readxl)

# --- data preparation (from otp_combined.R) ---
df_therapeutic_meta <- read_excel("metadata/otp.xlsx", sheet = "therapeutic_area")

df_disease_annot <- fread("results/otp/26.03/all_disease.tsv") %>% 
  left_join(df_therapeutic_meta, by = join_by(therapeuticAreaName == therapeutic_area)) %>% 
  filter(include) %>% 
  group_by(id, name) %>% 
  summarise(area_label = list(area_label))

df_chembl <- fread("results/otp/26.03/all_drug_evidence_chembl.tsv") %>% 
  .[score >= 0.1, .SD[which.max(score)], by = .(targetId, diseaseId)]

phases <- c("Phase I", "Phase II", "Phase III",
            "Phase I+", "Phase II+", "Phase III+", "Approved")

gene_universe_otp <- intersect(df_msmr_tenk10k[, unique(probeID)], df_chembl[, unique(targetId)])

disease_universe <- intersect(
  df_trait_map[supercategory == "disease", query_id],
  df_chembl[score >= 0.1, unique(diseaseId)]
)

df_msmr_tenk10k[df_trait_map, query_id := i.query_id, on = c("phenotype" = "trait_id")]

df_mr_max <- df_msmr_tenk10k %>%
  filter(probeID %in% gene_universe_otp & query_id %in% disease_universe) |> 
  group_by(probeID, query_id) |> 
  summarise(mr = max(mr, na.rm = TRUE),
            min_lfdr = min(lfdr_msmr_pheno, na.rm = TRUE),
            dir = sign(sum(sign(b_SMR) * mr, na.rm = TRUE)),
            min_p = min(p_SMR_multi))
setDT(df_mr_max)

df_mr_phase <- df_chembl %>% 
  filter(targetId %in% gene_universe_otp & 
           diseaseId %in% disease_universe) %>%
  group_by(targetId, diseaseId) %>% 
  slice_max(score) %>% 
  filter(score >= 0.1) %>%
  mutate(phase = case_when(
    score >= 0.1 & score < 0.2 ~ "Phase I",
    score >= 0.2 & score < 0.7 ~ "Phase II",
    score >= 0.7 & score < 0.8 ~ "Phase III",
    score >= 0.8 ~ "Approved"
  ) %>% factor(phases)
) %>% 
  setDT()

df_otp_assoc <- read_parquet("results/otp/26.03/otp_assoc_overall.gz.parquet") %>%
  filter(targetId %in% gene_universe_otp, diseaseId %in% disease_universe) %>%
  select(diseaseId, targetId, evidenceCount, association_type, associationScore) |> 
  pivot_wider(names_from = association_type, values_from = associationScore)
setDT(df_otp_assoc)

df_mr_phase[df_otp_assoc, `:=`(otp_assoc_direct = i.direct,
                               otp_assoc_indirect = i.indirect), 
            on = c("targetId", "diseaseId")]
df_mr_phase[, otp_assoc_max := pmax(otp_assoc_direct, otp_assoc_indirect, na.rm = TRUE)]

setDT(df_disease_annot)
df_mr_phase[df_disease_annot, `:=`(disease_label = i.name), on = c("diseaseId" = "id")]
df_mr_phase[df_mr_max, `:=`(mr = i.mr, min_p = i.min_p, direction = i.dir),
            on = c("targetId" = "probeID", "diseaseId" = "query_id")]
df_mr_phase[df_gene_annot, `:=`(hgnc_symbol = i.hgnc_symbol), on = c("targetId" = "ensembl_gene_id")]

# --- plot ---
# built at final print size for the journal (max 180 x 225 mm, text 5-7 pt)
fig_w <- 180 / 25.4
fig_h <- 225 / 25.4
sz <- list(label = 5, strip = 6.5, base = 7, point = 1.2, ncol = 5)

(p_top_drug_all <- df_mr_phase %>%
    filter(mr == 1) %>% 
    mutate(pheno_label = str_to_title(disease_label),
           score = otp_assoc_max) %>% 
    group_by(pheno_label) %>%
    mutate(rank_score = frank(-score)) %>% 
    ggplot(aes(x = -rank_score, y = score)) +
    theme_classic(base_size = sz$base) +
    geom_point(aes(fill = phase), shape = "circle filled", size = sz$point, stroke = 0.3) +
    facet_wrap(~pheno_label, scale = "free_x",
               ncol = sz$ncol, axes = "all",
               labeller = labeller(pheno_label = label_wrap_gen(20))) +
    theme(legend.position = "none") +
    geom_text_repel(
      aes(label = hgnc_symbol),
      fontface = "italic",
      size = sz$label/.pt, segment.size = 0.1,
      segment.color = "black",
      point.padding = 0.1,
      min.segment.length = 0,
      data = ~filter(.x, rank_score <= 15),
      max.overlaps = Inf,
      direction = "y",
      xlim = c(Inf, NA),
      force = 6,
      box.padding = 0.05,
      max.iter = 1e5,
      hjust = 1,
      seed = 42
    ) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0)),
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       guide = guide_axis(cap = "both")) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(y = "Association Score", x = NULL) +
    scale_fill_viridis_d(
      direction = -1,
      name = "Maximum clinical development stage",
      guide = guide_legend(override.aes = list(size = 2.5))) +
    theme(strip.clip = "off",
          strip.background = element_blank(),
          plot.margin = margin(r = 30, t = 2, l = 1, b = 1, unit = "pt"),
          panel.spacing.x = unit(32, "pt"),
          panel.spacing.y = unit(6, "pt"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.text = element_text(face = "bold", size = sz$strip,
                                    hjust = 0),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray90"),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          # legend in the empty space beside the last (single-panel) row
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.06),
          legend.justification.inside = c(1, 0.5),
          legend.direction = "horizontal",
          legend.title.position = "top")
)

# Extended Data Figure 8
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig-e8.pdf",
       p_top_drug_all, device = cairo_pdf, width = fig_w, height = fig_h, bg = "white")
ggsave("figures/fig-e8.png",
       p_top_drug_all, device = agg_png, width = fig_w, height = fig_h, dpi = 300, bg = "white")
