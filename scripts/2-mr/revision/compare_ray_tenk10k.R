# Compare eQTL MR results: Ray et al (OneK1K) vs TenK10K for CAD
library(paletteer)
library(data.table)
library(tidyverse)
library(qvalue)
library(ggrepel)
library(patchwork)
library(ragg)
library(scales)

source("scripts/preprocess_strict.R")

# ── 1. TenK10K: filter for CAD ──────────────────────────────────────────────

df_cad_tenk10k <- df_msmr_tenk10k[phenotype == "cad"]

# ── 2. Ray et al: load harmonised data and compute MR per gene × cell type ──

df_mr_ray <- fread("resources/misc/onek1k_cad_mr_ray2025.txt")
setnames(df_mr_ray, "Cell type", "cell_type_ray")

# IVW (multiple instruments) or Wald ratio (single instrument)
# compute_mr <- function(dt) {
#   n <- nrow(dt)
#   if (n == 0) return(NULL)
#   if (n == 1) {
#     b  <- dt$beta.outcome / dt$beta.exposure
#     se <- abs(dt$se.outcome / dt$beta.exposure)
#   } else {
#     # IVW: weight = (beta.exposure / se.outcome)^2
#     w  <- (dt$beta.exposure / dt$se.outcome)^2
#     b  <- sum(w * dt$beta.outcome / dt$beta.exposure) / sum(w)
#     se <- sqrt(1 / sum(w))
#   }
#   pval <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
#   list(b_ray = b, se_ray = se, p_ray = pval, nsnp_ray = n,
#        method_ray = if (n == 1) "Wald ratio" else "IVW")
# }

# mr_ray <- df_mr_ray[,
#   compute_mr(.SD),
#   by = .(Gene = exposure, cell_type_ray = cell_type)
# ]

# FDR correction across all gene × cell-type tests
df_mr_ray[, sig_ray  := qval < 0.05]

# ── 3. Cell type mapping: Ray et al → TenK10K major_cell_type ────────────────

# BM (bone marrow) has no direct equivalent in TenK10K PBMC data; excluded
cell_type_map <- tribble(
  ~cell_type_ray, ~major_cell_type,
  "B",     "B",
  "CD4ET", "CD4 T",
  "CD4NC", "CD4 T",
  "CD4S",  "CD4 T",
  "CD8ET", "CD8 T",
  "CD8NC", "CD8 T",
  "CD8S",  "CD8 T",
  "DC",    "Dendritic",
  "MonoC", "Monocyte",
  "MonoNC","Monocyte",
  "NK",    "NK",
  "NKR",   "NK",
  "Plasma","Plasma B"
) |> setDT()

df_mr_ray[cell_type_map, major_cell_type := i.major_cell_type, on = "cell_type_ray"]

# Aggregate Ray results to major cell type (take most significant per gene × major CT)
mr_ray_major <- df_mr_ray[
  !is.na(major_cell_type),
  .SD[which.min(qval)],
  by = .(Gene_ID, major_cell_type)
][, .(Gene = Gene_ID, major_cell_type, b_ray = b, se_ray = se, p_ray = pval, qvalue_ray = qval, nsnp_ray = nsnp, sig_ray)]

# ── 4. Aggregate TenK10K to major cell type ──────────────────────────────────

df_cad_major <- df_cad_tenk10k[
  !is.na(qval_msmr_pheno),
  .SD[which.min(qval_msmr_pheno)],
  by = .(Gene, major_cell_type)
][, .(Gene, probeID, major_cell_type, b_SMR, se_SMR, p_SMR,
      qval_msmr_pheno, p_HEIDI)]

df_cad_major[, sig_tenk10k := qval_msmr_pheno < 0.05]

# ── 5. Merge and classify overlap ────────────────────────────────────────────

df_compare <- merge(df_cad_major, mr_ray_major,
  by = c("Gene", "major_cell_type"), all = TRUE)

df_compare[is.na(sig_tenk10k), sig_tenk10k := FALSE]
df_compare[is.na(sig_ray),     sig_ray     := FALSE]

df_compare[, overlap := fcase(
  sig_tenk10k & sig_ray,  "TenK10K & Ray et al",
  sig_tenk10k & !sig_ray, "TenK10K only",
  !sig_tenk10k & sig_ray, "Ray et al only",
  default = "Neither"
)]

# ── 6. Summary counts ────────────────────────────────────────────────────────

cat("=== Significant gene × major-cell-type pairs ===\n")
cat("TenK10K (CAD):", sum(df_cad_major$sig_tenk10k, na.rm = TRUE), "\n")
cat("Ray et al:    ", sum(mr_ray_major$sig_ray, na.rm = TRUE), "\n")
df_overlap_n <- df_compare[sig_tenk10k | sig_ray, .N, by = overlap]
print(df_overlap_n)

# ── 7. Plot: overlap bar chart by major cell type ────────────────────────────

major_ct_order <- c("CD4 T", "CD8 T", "Unconventional T", "NK",
                    "Plasma B", "B", "Monocyte", "Dendritic")

df_plot_overlap <- df_compare[
  sig_tenk10k | sig_ray,
  .N,
  by = .(overlap, major_cell_type)
][major_cell_type %in% major_ct_order]

df_plot_overlap[, major_cell_type := factor(major_cell_type, levels = major_ct_order)]

overlap_colours <- c(
  "TenK10K & Ray et al"          = "#2ca25f",
  "TenK10K only"  = "#2171b5",
  "Ray et al only"= "#d95f02"
)

p_overlap <- ggplot(df_plot_overlap,
    aes(x = major_cell_type, y = N, fill = overlap)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = overlap_colours, name = NULL) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  labs(x = NULL, y = "Gene × cell-type pairs",
       title = "Overlap of significant MR findings for CAD")

# ── 8. Plot: beta concordance scatter (genes in Both) ────────────────────────

df_both <- df_compare[overlap == "TenK10K & Ray et al"]

# concordance statistics
n_concordant <- df_both[sign(b_SMR) == sign(b_ray), .N]
n_total_both <- nrow(df_both)

p_scatter <- ggplot(df_both, aes(x = b_ray, y = b_SMR)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "firebrick", linewidth = 0.6) +
  geom_point(aes(colour = major_cell_type), alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = Gene), size = 2.2, max.overlaps = 20, segment.size = 0.1, segment.color = "grey50", fontface = "italic") +
  scale_colour_paletteer_d("ggthemes::Tableau_10", name = "Major cell type") +
  theme_bw() +
  labs(
    x = bquote(beta[MR] ~ "— Ray et al (OneK1K)"),
    y = bquote(beta[MR] ~ "— TenK10K")
  )

# ── 9. Plot: -log10 p comparison (all shared genes) ─────────────────────────
df_shared_genes <- df_compare[
  Gene %in% df_compare[sig_tenk10k | sig_ray, unique(Gene)]
][!is.na(p_SMR) & !is.na(p_ray)]

p_pval <- ggplot(df_shared_genes[overlap != "Neither"],
    aes(x = -log10(p_ray), y = -log10(p_SMR), colour = overlap)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_label_repel(
    data = df_shared_genes[overlap == "TenK10K & Ray et al"],
    aes(label = Gene), size = 2, max.overlaps = Inf,
    box.padding = 0.3, point.padding = 0.2, segment.color = "grey50",
    fill = alpha("white", 0.7), label.size = 0, label.padding = 0.1,
    segment.size = 0.1, fontface = "italic", show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = overlap_colours, name = "MR association with CAD") +
  theme_bw() +
  labs(
    x = expression(-log[10]~italic(P) ~ "Ray et al (OneK1K)"),
    y = expression(-log[10]~italic(P) ~ "TenK10K")
  )

# ── 10. Save outputs ──────────────────────────────────────────────────────────

dir.create("figures/revision", showWarnings = FALSE, recursive = TRUE)
dir.create("results/revision",  showWarnings = FALSE, recursive = TRUE)

ggsave("figures/revision/compare_ray_tenk10k_overlap.png", p_overlap,
  device = agg_png, width = 7, height = 5, scaling = 0.8)

ggsave("figures/revision/compare_ray_tenk10k_scatter.png", p_scatter,
  device = agg_png, width = 8, height = 6, scaling = 1.2)

ggsave("figures/revision/compare_ray_tenk10k_pval.png", p_pval,
  device = agg_png, width = 8, height = 6, scaling = 1.2)

fwrite(df_compare[sig_tenk10k | sig_ray],
  "results/revision/compare_ray_tenk10k_significant.tsv", sep = "\t")
