# OTP overall evidence
source("scripts/preprocess_strict.R")

library(ggridges)
library(ggrepel)
library(scales)
library(patchwork)
library(ragg)
library(paletteer)
library(readxl)
library(broom)
library(sciscales)

df_therapeutic_meta <- read_excel("resources/metadata/otp.xlsx", sheet = "therapeutic_area")

df_phase_meta <- read_excel("resources/metadata/otp.xlsx", sheet = "phase") %>% 
  mutate(phase_label = factor(phase_label, unique(phase_label), ordered = TRUE))

df_disease_annot <- fread("results/otp/26.03/all_disease.tsv") %>% 
  left_join(df_therapeutic_meta, by = join_by(therapeuticAreaName == therapeutic_area)) %>% 
  filter(include) %>% 
  group_by(id, name) %>% 
  summarise(area_label = list(area_label))

# df_drug <- read_parquet("results/otp/26.03/otp_drug.gz.parquet")
# df_drug_summary <- read_parquet("results/otp/25.03/otp_drug_summary.gz.parquet")

df_mechanism <- fread("results/otp/26.03/all_drug_mechanism_chembl.tsv")
df_mechanism_info <- fread("resources/misc/chembl_action_type.tsv")
df_mechanism[df_mechanism_info, parent_mechanism := str_to_sentence(i.PARENT_TYPE),
             on = c("actionType" = "ACTION_TYPE")]

df_chembl <- fread("results/otp/26.03/all_drug_evidence_chembl.tsv") %>% 
  .[score >= 0.05, .SD[which.max(score)], by = .(targetId, diseaseId)]

df_chembl[df_mechanism, `:=`(mechanism = i.parent_mechanism), 
          on = c("targetId", "drugId")]

df_phase <- list(
  `New Drug` = expression(score >= 0.05 & score < 0.1),
  `Phase I` = expression(score >= 0.1 & score < 0.2),
  `Phase II` = expression(score >= 0.2 & score < 0.7),
  `Phase III` = expression(score >= 0.7 & score < 0.8),
  `Phase I+` = expression(score >= 0.1),
  `Phase II+` = expression(score >= 0.2),
  `Phase III+` = expression(score >= 0.7),
  `Approved` = expression(score >= 0.8)
) %>% 
  map_df(~df_chembl[eval(.x)], .id = "name") %>% 
  select(-score) %>% 
  mutate(value = TRUE) %>% 
  pivot_wider() %>% 
  setDT()

# gene_universe_otp <- df_msmr_tenk10k[gene_type == "protein_coding", unique(probeID)]
# gene_universe_otp <- df_chembl[, unique(targetId)]
gene_universe_otp <- intersect(df_msmr_tenk10k[,unique(probeID)], df_chembl[, unique(targetId)])

disease_universe <- intersect(
  df_trait_map[supercategory == "disease", query_id],
  df_phase[, unique(diseaseId)]
)

phases <- c("Phase I", "Phase II", "Phase III",
            "Phase I+", "Phase II+", "Phase III+", "Approved")
phases_mini <- c("Phase I", "Phase II", "Phase III", "Approved")

df_msmr_tenk10k[df_trait_map, query_id := i.query_id, on = c("phenotype" = "trait_id")]

df_mr_max <- df_msmr_tenk10k %>%
  filter(probeID %in% gene_universe_otp & query_id %in% disease_universe) |> 
  group_by(probeID, query_id) |> 
  # filter(mr == TRUE) |>
  summarise(mr = max(mr, na.rm = TRUE),
            min_lfdr = min(lfdr_msmr_pheno, na.rm = TRUE),
            dir = sign(sum(sign(b_SMR) * mr, na.rm = TRUE)),
            min_p = min(p_SMR_multi))
setDT(df_mr_max)
# fwrite(df_mr_max, "resources/misc/df_mr_max.tsv", sep = "\t")
# df_mr_max <- fread("resources/misc/df_mr_max.tsv")

gene_disease_universe <- expand_grid(gene = gene_universe_otp, disease = disease_universe)
# gene_disease_universe <- unique(df_phase[,.(gene = targetId, disease = diseaseId)])

df_gd <- gene_disease_universe %>% 
  left_join(df_phase, by = c(gene = "targetId", disease = "diseaseId")) %>% 
  mutate(across(c(all_of(phases)),
                ~nafill(as.numeric(.x), fill = FALSE) %>% as.logical))
setDT(df_gd)
df_gd[, mr := FALSE]
df_gd[df_mr_max, mr := i.mr, on = c(gene = "probeID", disease = "query_id")]
df_gd[df_mr_max, min_lfdr := i.min_lfdr, on = c(gene = "probeID", disease = "query_id")]

# df_gd <- df_phase |> 
#   left_join(df_mr_max, by = c(targetId = "probeID", diseaseId = "query_id")) |> 
#   mutate(across(c(all_of(phases), mr),
#                 ~nafill(as.numeric(.x), fill = FALSE) %>% as.logical))

# get numbers
df_gd[, .N, by = mr] %>% mutate(prop = N / sum(N))

# number of t-i pair available in chembl
df_gd_analysis <- df_gd[!is.na(drugId)]
df_gd_analysis[, .N]
df_gd_analysis[, n_distinct(gene)]
df_gd_analysis[, n_distinct(disease)]

# number of t-i pairs with MR support & percentage
# df_gd_analysis[mr == TRUE, .N]
# df_gd_analysis[, sum(mr == TRUE) / .N]

# plot
# df_gd_analysis |> 
#   pivot_longer(all_of(phases_mini), names_to = "phase", values_to = "value") |>
#   mutate(phase = factor(phase, levels = phases_mini)) |>
#   filter(value == TRUE) |>
#   ggplot(aes(y = -log10(min_p), x = phase)) +
#   geom_point(position = position_jitter(width = 0.2, height = 0)) +
#   theme_bw()

# Drug target enrichment analysis over a range of LFDR thresholds
# calc_enrichment <- function(lfdr_threshold){
#   df <- copy(df_gd)
#   df[is.na(min_lfdr), min_lfdr := 1]
#   df[, mr := min_lfdr < lfdr_threshold]
#   df_results <- map_df(phases, ~{
#    table(df$mr, df[[.x]]) |> fisher.test(alternative = "greater") |> tidy() |> 
#      mutate(phase = .x, lfdr_threshold = lfdr_threshold)
#   })
#   return(df_results)
# }

# thresholds <- c(0.005, 0.01, 0.05, 0.1, 0.2)

# df_threshold_enrichment <- map_df(thresholds, ~calc_enrichment(.x))


df_msmr_stats <- df_gd %>%
  pivot_longer(all_of(phases), names_to = "phase", values_to = "value") %>% 
  mutate(across(c(mr, value), ~factor(.x, levels = c("TRUE", "FALSE"))))

df_stats_overall <- df_msmr_stats %>% 
  select(gene, disease, mr, value, phase) %>% 
  group_by(phase) %>%
  nest(data = c(gene, disease, mr, value)) %>% 
  mutate(table = map(data, ~table(mr = .x$mr, phase = .x$value)),
         stats = map(table, ~fisher.test(.x, alternative = "greater") %>% tidy())) %>% 
  mutate(disease_label = "Overall",
         category = "Overall")

df_stats_overall %>% unnest(stats)

# per disease analysis
df_stats_disease <- df_gd %>%
  select(gene, disease, all_of(phases), mr) %>% 
  pivot_longer(all_of(phases), 
               names_to = "phase", values_to = "value") %>% 
  mutate(across(c(mr, value), ~factor(.x, levels = c("TRUE", "FALSE")))) %>% 
  left_join(df_disease_annot %>%
              select(disease = id, disease_label = name) %>% 
              distinct()) %>% 
  group_by(phase, disease, disease_label) %>%
  nest(data = c(gene, mr, value)) %>% 
  mutate(table = map(data, ~table(mr = .x$mr, phase = .x$value)),
         stats = map(table, ~fisher.test(.x, alternative = "greater") %>% tidy()))

df_all <- bind_rows(
  list(Overall = df_stats_overall, Disease = df_stats_disease),
  .id = "category"
) %>% 
  # filter min 5 obs in each cell within table
  # filter(map_lgl(table, ~all(.x > 5))) %>% 
  select(phase, category, disease, disease_label, stats) %>% 
  unnest(stats) %>% 
  filter(p.value < 1)
  # group_by(phase) %>% 
  # mutate(qvalue = qvalue_truncp(p.value)$qvalues)

# bar chart overall
tally_gene <- df_chembl %>% 
  filter(targetId %in% gene_universe_otp & 
         diseaseId %in% disease_universe) %>%
  filter(score >= 0.1) %>%
  group_by(targetId) %>% 
  summarise(score = max(score)) %>% 
  mutate(phase = case_when(
    score == 0.1 ~ "Phase I",
    score == 0.2 ~ "Phase II",
    score == 0.7 ~ "Phase III",
    TRUE ~ "Approved"
  ) %>% factor(phases),
  mr_genes = targetId %in% df_mr_max[mr == 1, unique(probeID)] %>% 
    factor(levels = c("TRUE", "FALSE"),
           labels = c("MR support", "No MR support"))
)

pals <- paletteer_d("MexBrewer::Frida")
(p_bar_overall <- tally_gene %>% 
  group_by(phase, mr_genes) %>% 
  tally() %>% 
  arrange(phase, mr_genes) %>% 
  ggplot(aes(x = phase, fill = mr_genes, y = n, group = phase)) +
  theme_bw() +
  geom_col(position = "stack") +
  scale_fill_manual(values = pals[3:2],
                    breaks = c("No MR support", "MR support")) +
  labs(x = "Maximum clinical development phase",
       y = "Number of targets",
       fill = NULL) +
  geom_label(aes(label = n), 
             position = position_stack(vjust = 0.5), 
             size = 9/.pt, color = "black", fill = alpha("white", 0.4),
             linewidth = 0,
             show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.line = element_line(),
        panel.border = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        legend.key.size = unit(0.75, "lines"),
        legend.key.spacing.y = unit(0.25, "lines"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
  )

# forest plot overall
get_sensitivity <- function(tb){
  tb[1,1] / sum(tb[,1])
}
get_sensitivity_text <- function(tb){
  paste(tb[1,1], "/", sum(tb[,1]))
}

(
p_enrich <- df_stats_overall %>% 
  # filter(phase %in% c("Phase I+", "Phase II+", "Phase III+", "Approved")) %>% 
  unnest(stats) %>% 
  mutate(phase = factor(phase, levels = rev(phases)),
         sens_text = map_chr(table, get_sensitivity_text)) %>% 
  ggplot(aes(y =  phase, x = estimate)) +
  theme_bw() +
  geom_point(shape = "diamond", size = 3, color = pals[4]) +
  geom_segment(aes(x = conf.low, xend = estimate, yend = phase),
               linewidth = 1, color = pals[4]) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Odds ratio", y = NULL) +
  scale_x_continuous(trans = "log", breaks = c(1, 1.2, 1.5, 2, 3),
                     guide = guide_axis(cap = "upper"), 
                     expand = expansion(mult = c(0.05, 0.6))) +
  geom_text(aes(label = number(estimate, .1)),
            hjust = 0.5, vjust = -1, size= 8 / .pt) +
  geom_text(aes(label = sens_text, x = Inf), 
            hjust = 1, size = 8 / .pt, color = "black") +
  annotate("text", y = Inf, x=Inf, vjust = 1, hjust = 1,
           color = "black", size = 8/.pt, fontface = "bold",
           label = "MR support / Total") +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = "plain"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r =1, unit = "lines"),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.major.y = element_blank())
)

# per disease results
df_mr_phase <- df_chembl %>% 
  filter(targetId %in% gene_universe_otp & 
           diseaseId %in% disease_universe) %>%
  group_by(targetId, diseaseId) %>% 
  slice_max(score) %>% 
  filter(score >= 0.1) %>%
  mutate(phase = case_when(
    score == 0.1 ~ "Phase I",
    score == 0.2 ~ "Phase II",
    score == 0.7 ~ "Phase III",
    TRUE ~ "Approved"
  ) %>% factor(phases)
) %>% 
  setDT()

df_otp_assoc <- read_parquet("results/otp/26.03/otp_assoc_overall.gz.parquet") %>%
  filter(targetId %in% gene_universe_otp, diseaseId %in% disease_universe) %>%
  pivot_wider(names_from = association_type, values_from = associationScore)
setDT(df_otp_assoc)

df_mr_phase[df_otp_assoc, `:=`(otp_assoc_direct = i.direct,
                               otp_evidence_count = i.evidenceCount,
                               otp_assoc_indirect = i.indirect), 
            on = c("targetId", "diseaseId")]
df_mr_phase[, otp_assoc_max := pmax(otp_assoc_direct, otp_assoc_indirect, na.rm = TRUE)]

setDT(df_disease_annot)
df_mr_phase[df_disease_annot,  `:=`(disease_label = i.name), 
            on = c("diseaseId" = "id")]
df_mr_phase[df_mr_max, `:=`(mr = i.mr, min_p = i.min_p, direction = i.dir),
            on = c("targetId" = "probeID", "diseaseId" = "query_id")]
df_mr_phase[df_gene_annot, `:=`(hgnc_symbol = i.hgnc_symbol),  on = c("targetId" = "ensembl_gene_id")]

df_mr_phase[, `:=`(mr_any = targetId %in% df_mr_max[mr == 1, unique(probeID)])]

# Open targets by datatype
gene_universe_tenk10k <- df_msmr_tenk10k[, unique(probeID)]
disease_universe_tenk10k <- df_trait_map[supercategory == "disease", unique(query_id)]

# genetic sources
# see: https://platform-docs.opentargets.org/associations
otp_data_sources <- list(
  genetic = c("gwas_credible_sets", "gene_burden", "clingen", "eva", "uniprot_literature", "uniprot_variants",
"gene2phenotype", "genomics_england", "orphanet"),
  non_genetic = c("europepmc", "clinical_precedence", "cancer_gene_census", "intogen", "eva_somatic", "cancer_biomarkers", "crispr_screen", "project_score", "reactome", "expression atlas", "impc")
)


df_otp_datasource <- read_parquet("results/otp/26.03/otp_assoc_datasource.gz.parquet") %>% 
  filter(targetId %in% gene_universe_otp, diseaseId %in% disease_universe) %>% 
  mutate(data_source = ifelse(aggregationValue %in% otp_data_sources$genetic, "genetic", "non_genetic")) %>% 
  pivot_wider(names_from = "association_type", values_from = "associationScore") %>% 
  mutate(score = pmax(direct, indirect, na.rm = TRUE)) %>% 
  group_by(data_source, targetId, diseaseId) %>% 
  summarise(score = mean(score, na.rm = TRUE)) %>% 
  pivot_wider(names_from = data_source, values_from = score)
setDT(df_otp_datasource)

df_mr_phase[df_otp_datasource, `:=`(otp_genetic = i.genetic, otp_nongenetic = i.non_genetic), 
            on = c("targetId", "diseaseId")]

# plot otp assoc - mr results
df_ttest <- df_mr_phase %>%
  select(mr, otp_assoc_max, otp_genetic, otp_nongenetic) %>%
  pivot_longer(c(otp_assoc_max, otp_genetic, otp_nongenetic), names_to = "source") %>% 
  filter(!is.na(value), !is.na(mr)) %>% 
  group_by(mr, source) %>% 
  summarise(value = list(value)) %>% 
  pivot_wider(names_from = "mr") %>% 
  mutate(ttest = map2(`0`, `1`, ~tidy(t.test(.x, .y))),
         wilcox = map2(`0`, `1`, ~tidy(wilcox.test(.x, .y))),
         median_no_mr = map_dbl(`0`, ~median(.x, na.rm = TRUE)),
         median_mr = map_dbl(`1`, ~median(.x, na.rm = TRUE))) %>% 
  unnest(wilcox) %>% 
  mutate(plab = to_scientific(p.value) %>% as.character)

(p_violin_mr <- df_mr_phase %>% 
  filter(!is.na(mr)) %>% 
  mutate(x_lab = ifelse(mr == 0, "No MR support", "MR support") %>% 
           factor(levels = c("No MR support", "MR support"))) %>% 
  select(x_lab, otp_assoc_max, otp_genetic, otp_nongenetic) %>%
  pivot_longer(-x_lab, names_to = "source") %>% 
  ggplot(aes(x = x_lab, y = value)) +
  theme_classic() +
  geom_violin(aes(fill = x_lab), alpha = 0.4) +
  geom_boxplot(aes(fill = x_lab), outliers = FALSE, width = 0.1) +
  facet_grid(cols = vars(source), labeller = as_labeller(
    c(otp_assoc_max = "Overall Association Score",
      otp_genetic = "Genetic Association Score",
      otp_nongenetic = "Non-Genetic Association Score")
  )) +
  annotate("errorbar",
    xmin = 1, xmax = 2,  y = 1.05, width = 0.01,
    color = "gray30", linewidth = 0.5, linetype = "solid", orientation = "y") +
  geom_text(
    aes(label = paste("italic(P)['median diff.'] ==", plab)),
    data = df_ttest, size = 9/.pt, color = "gray30",
    y = 1.05, vjust = -0.5, x = 1.5, parse = TRUE) +
  scale_fill_manual(values = pals[c(2,3)], guide = "none", aesthetics = c("color", "fill")) +
  scale_y_continuous(expand = expansion(add = c(0.01, 0.12)),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     guide = guide_axis(cap = "both")) +
  labs(x = NULL, y = "Association Score") +
  coord_cartesian(clip = "off") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        strip.text  = element_text(face = "bold", size = 9),
        panel.spacing = unit(c(1), "lines"),
        panel.grid.major.y = element_line(color = "gray90"))
)

# Show drug by mechanism of action
pals2 <- paletteer_d("RColorBrewer::Spectral")

# calculate agreement between observed (ChEMBL) and predicted (MR) direction of effect
df_metrics <- df_mr_phase |> 
    filter(mr == 1, direction != 0, mechanism %in% c("Positive modulator", "Negative modulator")) |> 
    mutate(observed_dir = ifelse(mechanism == "Positive modulator", 1, -1),
           predicted_dir = -direction) |> 
    select(targetId, diseaseId, phase, observed_dir, predicted_dir)

# calculate P-value for concordance by phase
#permutation test 
set.seed(123)
n_perms <- 10000
permuted_scores <- numeric(n_perms)
for(i in 1:n_perms) {
  shuffled_pred <- sample(df_metrics$predicted_dir) # Shuffle
  permuted_scores[i] <- mean(df_metrics$observed_dir == shuffled_pred)
}

df_concordance <- df_metrics |> 
  group_by(phase) |>
  summarise(concordance = mean(observed_dir == predicted_dir),
            n_concordant = sum(observed_dir == predicted_dir),
            n_discordant = sum(observed_dir != predicted_dir),
            n = n()) |>
  group_by(phase) |>
  mutate(p_permutation = (sum(permuted_scores >= concordance)  + 1) / (n_perms + 1))

# visualise results
(p_bar_concordance <- df_concordance |> 
  pivot_longer(c(n_concordant, n_discordant), names_to = "type", values_to = "count") |>
  mutate(type = factor(type, levels = c("n_discordant", "n_concordant"))) |> 
  arrange(phase, type) |>
  mutate(prop = count / sum(count), x_label = cumsum(count) - (count / 2)) |>
  ggplot(aes(y = fct_rev(phase), x = count)) +
  geom_col(aes(fill = fct_rev(type)), position = "stack", color = "black", width = 0.75) +
  scale_fill_manual(values = pals[c(3,2)],
    labels = c(n_discordant ="MoA not matching MR direction of effect",
               n_concordant = "MoA matching MR direction of effect"),
    guide = guide_legend(reverse = TRUE),
    name = NULL,aesthetics = c("fill")) +
  # Add P-value on the right
  geom_text(
    aes(x = Inf, label = paste0(to_scientific(p_permutation))),
    data = df_concordance, 
    hjust = 1, vjust = 0.5, size = 9/.pt, color = "black", parse = TRUE
  ) +
  annotate("text", x = Inf, y = Inf, vjust = 1, hjust = 1,
       color = "black", size = 9/.pt,
       label = "bold(italic(P)[permutation])", parse = TRUE) +
  # add % for each category
  geom_label(aes(x = x_label, label = percent(prop, 1)),
              hjust = 0.5, vjust = 0.5, size = 9/.pt, color = "black",
              fill = alpha("white", 0.4), linewidth = 0
  ) +
  labs(x = "Count", y = NULL) +
  scale_x_continuous(limits = c(0, 40),
                     breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0.0, 0.3))) +
  guides(x = guide_axis(cap = "upper")) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.5, "lines"),
        legend.direction = "vertical",
        legend.title = element_text(size = 9),
        legend.text = element_text(margin = margin(l = 0.5, unit= "lines")),
        legend.justification.bottom = "left",
        legend.key.spacing.y = unit(0.25, "lines"),
        legend.box.margin = margin(t = -0.5, unit = "lines"),
        strip.clip = "off",
        strip.text  = element_text(face = "bold", size = 9),
        panel.spacing = unit(c(1), "lines"),
        panel.grid.major.x = element_line(color = "gray90"))
)

ggsave("figures/strict/otp/p_bar_concordance.png", p_bar_concordance, scaling = 1.5, device = agg_png,
       width = 8, height = 6)
# (
# p_direction <- df_mr_phase %>%
#   filter(mr == 1, direction != 0, mechanism %in% c("Positive modulator", "Negative modulator")) %>% 
#   mutate(matched_mr_dir = (direction == 1 & mechanism == "Negative modulator") | (direction == -1 & mechanism == "Positive modulator")) %>% 
#   group_by(mechanism, matched_mr_dir, phase) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   complete(matched_mr_dir, nesting(mechanism, phase), fill = list(n = 0)) %>% 
#   ggplot(aes(x = as.factor(mechanism), y = n)) +
#   geom_col(aes(fill = as.factor(matched_mr_dir)), position = "dodge", color = "black", width = 0.7) +
#   facet_wrap(~phase, nrow = 1, scale = "free_x", space = "free_x") +
#   scale_fill_manual(
#     values = c(`TRUE` = pals2[2], `FALSE` = rev(pals2)[2]),
#     labels = c(`TRUE` = "Matched direction", `FALSE` = "Unmatched direction"),
#     breaks = c("TRUE", "FALSE"),
#     name = "MR effect direction"
#   ) +
#   labs(x = NULL, y = "Count") +
#   theme_classic() +
#   scale_y_continuous(expand = expansion(mult = c(0,0.02))) +
#   theme(axis.text.x = element_text(size = 8, angle = 40, hjust = 1),
#         strip.background = element_blank(),
#         strip.clip = "off",
#         strip.text = element_text(face = "bold"),
#         legend.position = "inside",
#         axis.title = element_text(size = 9),
#         legend.position.inside = c(0.01,0.99),
#         legend.justification.inside = c(0,1),
#         legend.title = element_text(size = 9),
#         panel.grid.major.y = element_line(color = "gray90"),
#         legend.key.size = unit(0.5, "lines"))
#   )

# Show top drug by association score
pos <- position_jitter(width = 0.5, height=0, seed = 12)
(p_mr_phase <- df_mr_phase %>% 
    filter(mr == 1) %>% 
    group_by(phase) %>% 
    mutate(rank_p = frank(min_p)) %>% 
    ggplot(aes(x = -rank_p, y = -log10(min_p))) +
    theme_classic() +
    geom_point(aes(fill = phase), size = 1.5, shape = "circle filled") +
    scale_fill_viridis_d(direction = -1, name = "Maximum clinical development stage") +
    facet_wrap(~phase, scale = "free_x", space = "free_x", nrow = 1) +
    coord_cartesian(clip = "off") +
    labs(y = bquote(-log[10]~italic(P)[MR]), x = NULL) +
    theme(strip.background = element_blank(),
          strip.clip = "off",
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 9),
          panel.grid.major.y = element_line(color = "gray90"),
          legend.position = "none",
          legend.background = element_rect(fill = NA),
          legend.key.spacing = unit(0.1, 'lines'),
          legend.position.inside = c(0,1),
          plot.margin = margin(l = 0.5, unit = "lines"),
          legend.justification.inside = c(0,1),
          strip.text = element_blank())
)

selected_diseases <- c("Crohn's disease" = "Crohn's disease", 
                       "rheumatoid arthritis" = "Rheumatoid arthritis",
                       "psoriasis" = "Psoriasis",
                       "type 2 diabetes mellitus" = "Type 2 diabetes",
                       "Alzheimer disease" = "Alzheimer's disease")
(p_top_drug <- df_mr_phase %>%
    filter(mr == 1) %>% 
    mutate(pheno_label = recode(disease_label, !!!selected_diseases),
           score = otp_assoc_max) %>% 
    filter(disease_label %in% names(selected_diseases)) %>%
    group_by(pheno_label) %>%
    mutate(rank_score = frank(-score)) %>% 
    ggplot(aes(x = -rank_score, y = score)) +
    theme_classic() +
    geom_point(aes(fill = phase), shape = "circle filled") +
    facet_wrap(~pheno_label, scale = "free_x",
               nrow = 1,  axes = "all",
               labeller = labeller(pheno_label = label_wrap_gen(25))) +
    theme(legend.position = "none") +
    geom_text_repel(
      aes(label = hgnc_symbol),
      fontface = "italic",
      size = 7/.pt, segment.size = 0.2,
      point.padding = 0.1,
      min.segment.length = 0,
      data = ~filter(.x, rank_score <= 5),
      max.overlaps = Inf,
      direction = "y",
      xlim = c(Inf, NA),
      force = 5,
      hjust = 1,
      segment.color = "black",
    ) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(add = c(0, 0.03)),
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       guide = guide_axis(cap = "both")) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(y = "Association Score", x = NULL) +
    scale_fill_viridis_d(
      direction = -1,
      name = "Maximum clinical development stage",
      guide = guide_legend(override.aes = list(size = 2))) +
    theme(strip.clip = "off",
          strip.background = element_blank(),
          plot.margin = margin(r = 2.5, t = 0, l = 0.2, unit = "lines"),
          panel.spacing.x = unit(3, "lines"),
          strip.text = element_text(face = "bold", size = 8,
                                    hjust = 0),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray90"),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom")
)

# combine
p_design <- "
AAABBB
AAABBB
CCCCCC
CCCCCC
DDDEEE
DDDEEE
DDDEEE
FFFFFF
FFFFFF
"

inside_tag <- theme(
  plot.tag.position = c(0,1),
  plot.tag = element_text(face = "bold", size = 14, hjust = 0, vjust = 1))

plots <- list(
  A = p_bar_overall + theme(legend.box.margin = margin(l = -1, unit= "lines")),
  B = wrap_elements(full = p_enrich),
  C = p_violin_mr,
  D = wrap_elements(full = p_bar_concordance),
  E = wrap_elements(plot = p_mr_phase),
  F = p_top_drug
) %>% 
  wrap_plots(guides = "keep", design = p_design)

plots <-  (plots &
  theme(
    plot.tag.position = c(0,1),
    plot.tag = element_text(face = "bold", size = 13, hjust = 0, vjust = 1)
    )) +
  plot_annotation(tag_levels = "a") &
  theme(axis.title = element_text(size = 9))


ggsave("figures/strict/otp/otp_combined_v5.png", device = agg_png, bg = "white",
       plots, width = 9, height = 12, dpi = 300, scaling = 1.2)

ggsave("figures/publication_pdf/Fig4 - Drug Identification Support.pdf",
       device = cairo_pdf, bg = "white",
       plots, width = 9, height = 12, dpi = 320, scale = 1/1.2)


# stand alone selected plot
ggsave("figures/strict/otp/otp_enrichment.png", device = agg_png, bg = "white",
       p_enrich, width = 5, height = 3, dpi = 300, scaling = 1.2)
ggsave("figures/strict/otp/selected_targets.png", device = agg_png, bg = "white",
       p_top_drug, width = 8, height = 4, dpi = 300, scaling = 1.2)


(p_top_drug_all <- df_mr_phase %>%
    filter(mr == 1) %>% 
    mutate(pheno_label = str_to_title(disease_label),
           score = otp_assoc_max) %>% 
    group_by(pheno_label) %>%
    mutate(rank_score = frank(-score)) %>% 
    ggplot(aes(x = -rank_score, y = score)) +
    theme_classic() +
    geom_point(aes(fill = phase), shape = "circle filled") +
    facet_wrap(~pheno_label, scale = "free_x",
               ncol = 4, axes = "all",
               labeller = labeller(pheno_label = label_wrap_gen(25))) +
    theme(legend.position = "none") +
    geom_text_repel(
      aes(label = hgnc_symbol),
      fontface = "italic",
      size = 5/.pt, segment.size = 0.1,
      segment.color = "black",
      
      point.padding = 0.1,
      min.segment.length = 0,
      data = ~filter(.x, rank_score <= 15),
      max.overlaps = Inf,
      direction = "y",
      xlim = c(Inf, NA),
      force = 2,
      hjust = 1,
      segment.color = "black",
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
      guide = guide_legend(override.aes = list(size = 2))) +
    theme(strip.clip = "off",
          strip.background = element_blank(),
          plot.margin = margin(r = 2.5, t = 0, l = 0.2, unit = "lines"),
          panel.spacing.x = unit(3, "lines"),
          strip.text = element_text(face = "bold", size = 8,
                                    hjust = 0),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray90"),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom")
)


ggsave("figures/strict/otp/otp_ti_mr_all.png", device = agg_png, bg = "white",
       p_top_drug_all, width = 8, height = 12, dpi = 300, scaling = 1.1)

# write supplementary tables
source("scripts/util/helper.R")


# table of max MR evidence & max clinical dev per target
tbl_tally_gene <- tally_gene %>% 
  left_join(df_gene_annot[,.(targetId = ensembl_gene_id, hgnc_symbol)]) %>% 
  arrange(desc(phase), mr_genes)

write_gs(tbl_tally_gene, "target_otp", 11)

# table of TI support estimate
tbl_stats <- df_stats_overall %>% 
  mutate(Evidence_MR = map_dbl(table, ~.x[1,1]),
         noEvidence_MR = map(table, ~.x[1,2]),
         Evidence_noMR = map(table, ~.x[2,1]),
         noEvidence_noMR = map(table, ~.x[2,2])) %>% 
  unnest(stats)

write_gs(tbl_stats, "ti_support_stats", 12)

# table of otp score by source of evidence
label <- tibble(
  source = c("otp_assoc_max", "otp_genetic", "otp_nongenetic"),
  label = c("Overall Association Score", "Genetic Association Score", "Non-Genetic Association Score")
)
tbl_otp_by_evidence <- df_ttest %>% 
  mutate(n_not_supported = map_dbl(`0`, ~length(.x)),
         n_supported = map_dbl(`1`, ~length(.x))) %>% 
  left_join(label) %>% 
  rename(source_id = source, source = label)
write_gs(tbl_otp_by_evidence, "otp_by_evidence", 13)


# table of TI support with MR

df_drug_molecule_chembl <- fread("results/otp/26.03/all_drug_molecule_chembl.tsv") %>% 
  rename(drug_name = name)

tbl_ti <- df_mr_phase %>% 
  filter(mr == 1) %>% 
  left_join(df_drug_molecule_chembl) %>% 
  left_join(df_metrics) |> 
  # add direction of effect
  arrange(desc(phase)) %>% 
  mutate(dir_label = case_when(
    direction == 1 ~ "Positive",
    direction == 0 ~ "Undetermined",
    direction == -1 ~ "Negative"),
  concordance = case_when(
    observed_dir == predicted_dir ~ "Concordant",
    observed_dir != predicted_dir ~ "Discordant",
    TRUE ~ "Undetermined"
  ))
write_gs(tbl_ti, "ti_drug_mr", 14)

