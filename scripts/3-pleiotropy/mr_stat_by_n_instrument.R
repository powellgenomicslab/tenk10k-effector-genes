# distribution of MR test statistics, evidence-tier counts, and a "fair" statistical
# comparison of single- vs multi-instrument MR results, to address the concern that
# single-instrument SMR tests may be less valid than multi-instrument ones.
suppressPackageStartupMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(arrow)
  library(ragg)
  library(scales)
  library(patchwork)
  library(paletteer)
}))

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/revision", recursive = TRUE, showWarnings = FALSE)

# source() runs at top level, so df_msmr_tenk10k AND df_msmr_eqtlgen (the eQTLGen
# MSMR results used for the eqtlgen_mr replication flag) both end up in this
# environment - the latter is reused below to build a fair replication test.
source("scripts/0-preprocess/preprocess_results.R")

setDT(df_msmr_tenk10k)

df_msmr_tenk10k[, instrument_cat := fifelse(n_snps_smr == 1, "Single instrument", "Multiple instruments (>1)")]
df_msmr_tenk10k[, instrument_cat := factor(instrument_cat,
  levels = c("Single instrument", "Multiple instruments (>1)")
)]

instrument_color <- paletteer_d("ggthemes::Tableau_10")[1:2] |>
  set_names(c("Single instrument", "Multiple instruments (>1)"))

# how many tests fall into each instrument category
df_n_tested <- df_msmr_tenk10k[, .(n_test = .N), by = instrument_cat]

# ---------------------------------------------------------------------------
# Panel a: distribution of the number of instruments (n_snps_smr), all tested
# pairs vs MR-significant hits only, to see whether the significant subset is
# enriched/depleted for any particular instrument count relative to the full
# tested population.
# ---------------------------------------------------------------------------
instrument_n_breaks <- c(0, 1, 2, 3, 4, 9, Inf)
instrument_n_labels <- c("1", "2", "3", "4", "5-9", "10+")
df_msmr_tenk10k[, instrument_n_bin := cut(n_snps_smr,
  breaks = instrument_n_breaks, labels = instrument_n_labels
)]

df_instrument_dist <- rbindlist(list(
  df_msmr_tenk10k[, .(n = .N, population = "All tests"), by = instrument_n_bin],
  df_msmr_tenk10k[mr == TRUE, .(n = .N, population = "MR-significant hits"), by = instrument_n_bin]
))
df_instrument_dist[, total := sum(n), by = population]
df_instrument_dist[, prop := n / total]
df_instrument_dist[, population := factor(population, levels = c("All tests", "MR-significant hits"))]

population_color <- c("All tests" = "#A4ABB0FF", "MR-significant hits" = "#4C6C94FF")

(p_panel_a <- ggplot(df_instrument_dist, aes(x = instrument_n_bin, y = prop, fill = population)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  geom_text(
    aes(label = percent(prop, accuracy = 0.1)),
    position = position_dodge(width = 0.8), vjust = -0.3, size = 8/.pt
  ) +
  scale_fill_manual(
    values = population_color,
    labels = c("All tests" = "All tests", "MR-significant hits" = bquote(LFDR[MR] < 0.05)),
    name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0.15)), labels = percent) +
  labs(
    x = "Number of instruments",
    y = "Proportion of tests",
    # title = "Distribution of instrument count",
    tag = "a"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "inside",
        legend.justification = c(1,1),
        panel.grid.major.y = element_line(color = "gray95"),
        legend.position.inside = c(0.99,0.99),
        legend.key.size = unit(0.75,"lines"))
)

# ---------------------------------------------------------------------------
# Panel b: evidence-tier counts (left) and proportions (right), side by side
#    (MR / MR+sensitivity / MR+sensitivity+coloc), stratified by instrument
#    category. Shown side by side rather than on a shared dual axis because the
#    count -> proportion conversion uses a DIFFERENT denominator per instrument
#    category (n_test), so a single secondary axis could not represent both
#    groups' proportions correctly at once.
#    Caveat: "sensitivity" = HEIDI pass OR Cochran's Q (IVW het.) pass OR MRLink2
#    pleiotropy pass. Cochran's Q (phet_ivw) needs >=2 instruments and is
#    therefore NA for 100% of single-instrument tests - single-instrument rows
#    can only pass "sensitivity" via HEIDI or MRLink2. This mechanically gives
#    multi-instrument tests one extra (OR-combined) chance to pass, so the raw
#    tier comparison below should NOT be read as evidence about validity on its
#    own - see the fair, criterion-matched comparison in section 4.
# ---------------------------------------------------------------------------
evidence_label <- c(
  mr = "MR",
  mr_sens = "MR + Sensitivity",
  mr_sens_coloc = "MR + Sensitivity + Coloc"
)

df_evidence_n <- df_msmr_tenk10k |>
  summarise(
    mr = sum(mr, na.rm = TRUE),
    mr_sens = sum(mr_sens, na.rm = TRUE),
    mr_sens_coloc = sum(mr_sens_coloc, na.rm = TRUE),
    n_test = n(),
    .by = instrument_cat
  ) |>
  pivot_longer(all_of(names(evidence_label)), names_to = "evidence", values_to = "n") |>
  mutate(
    evidence = factor(evidence, levels = names(evidence_label), labels = evidence_label),
    prop = n / n_test
  )

# statistical test: does the proportion of ALL tests reaching each evidence
# tier differ between single- and multi-instrument tests? Two-proportion
# z-test (prop.test) per tier, with the odds ratio as effect size - with N in
# the millions, p-values alone are not informative (all three will be tiny
# regardless of how big the underlying difference actually is).
# evidence_cols <- names(evidence_label)
# df_tier_test <- rbindlist(lapply(evidence_cols, function(col) {
#   tab <- df_msmr_tenk10k[, .(x = sum(get(col), na.rm = TRUE), n_test = .N), by = instrument_cat]
#   x_single <- tab[instrument_cat == "Single instrument", x]
#   n_single <- tab[instrument_cat == "Single instrument", n_test]
#   x_multi <- tab[instrument_cat == "Multiple instruments (>1)", x]
#   n_multi <- tab[instrument_cat == "Multiple instruments (>1)", n_test]
#   pt <- prop.test(c(x_multi, x_single), c(n_multi, n_single))
#   OR <- (x_multi / (n_multi - x_multi)) / (x_single / (n_single - x_single))
#   data.table(
#     evidence = col,
#     prop_single = x_single / n_single, prop_multi = x_multi / n_multi,
#     p_value = pt$p.value, OR = OR
#   )
# }))
# df_tier_test[, evidence := factor(evidence, levels = evidence_cols, labels = evidence_label)]
# # place all three annotations at one shared height (not per-tier max) - the
# # tiers span a >8-fold range in bar height, so a per-tier-relative offset puts
# # the smallest tier's annotation right on top of its own (tiny) bars/labels
# df_tier_test[, y_annot := max(df_evidence_n$prop) * 1.18]
# df_tier_test[, label := paste0(
#   "OR=", sprintf("%.2f", OR), ", ",
#   scales::pvalue(p_value, accuracy = 0.001, add_p = TRUE)
# )]

# print(df_tier_test)
# fwrite(df_tier_test, "results/revision/mr_instrument_evidence_tier_test.tsv", sep = "\t")

(
p_evidence_n <- ggplot(df_evidence_n, aes(x = evidence, y = n, fill = instrument_cat)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  geom_text(
    aes(label = comma(n)),
    position = position_dodge(width = 0.8), vjust = -0.3, size = 7/.pt
  ) +
  scale_fill_manual(values = instrument_color, name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0.15)), labels = comma) +
  scale_x_discrete(labels = ~str_wrap(.x, width = 16), expand = expansion(add =0.5)) +
  labs(x = NULL, y = NULL,
       title = "Number of significant associations",
      #  title = "Evidence-tier counts",
       tag = "b") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11),
  legend.title = element_text(size = 11),
          panel.grid.major.y = element_line(color = "gray95"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7.5))
)
(p_evidence_prop <- ggplot(df_evidence_n, aes(x = evidence, y = prop, fill = instrument_cat)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  geom_text(
    aes(label = percent(prop, accuracy = 0.1)),
    position = position_dodge(width = 0.8), vjust = -0.3, size = 8/.pt
  ) +
  # geom_text(
  #   data = df_tier_test, aes(x = evidence, y = y_annot, label = label),
  #   inherit.aes = FALSE, size = 6.5/.pt
  # ) +
  scale_fill_manual(values = instrument_color, name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0.4)), labels = percent) +
  scale_x_discrete(labels = ~str_wrap(.x, width = 16), expand = expansion(add = 0.5)) +
  labs(
    x = NULL, y = NULL,
    title = "Proportion of tests that are significant"
    # title = "Evidence-tier proportions",
    # caption = "Note: Cochran's Q (phet_ivw) requires at least 2 instruments and is NA for all\nsingle-instrument tests, so 'Sensitivity' gives multi-instrument tests an extra\nOR-combined pass criterion. See the fair, criterion-matched test in the next figure."
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "inside",
    legend.justification = c(1,1),
    legend.position.inside = c(0.99, 0.99),
    panel.grid.major.y = element_line(color = "gray95"),
    plot.title = element_text(size = 11),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7.5),
    # plot.caption = element_text(hjust = 0, size = 7, color = "grey30")
  )
)

# ---------------------------------------------------------------------------
# Panel c: QQ plot of the MR test p-value (p_SMR_multi) against its uniform null,
#    stratified by instrument category, over ALL tested pairs (not just
#    mr == TRUE - a QQ plot's whole point is to show departure from the null
#    across the full range, with the tail revealing excess significant hits).
#    If single- and multi-instrument MR are finding similar signal, the two
#    curves should track each other; if one departs from the diagonal harder/
#    earlier than the other, that group is yielding more (or stronger)
#    discoveries than the null would predict.
#    7.5M points can't be rendered directly: the most significant points are
#    kept in full and the bulk is log-spaced subsampled, which preserves the
#    curve's shape (a monotonic sorted sequence) without plotting every point.
#    Sits beside panel d; each keeps its own legend (plot_layout guides = "keep")
#    since the two encode different things (instrument category vs test type).
# ---------------------------------------------------------------------------
thin_qq_idx <- function(n, n_keep_top = 50000, n_bulk = 5000) {
  if (n <= n_keep_top) return(seq_len(n))
  top <- seq_len(n_keep_top)
  bulk <- unique(round(exp(seq(log(n_keep_top + 1), log(n), length.out = n_bulk))))
  sort(unique(c(top, bulk)))
}

make_qq_df <- function(p, group_label) {
  # clamp rather than drop exact-zero p-values (float underflow) - dropping
  # would discard exactly the most significant points, which matter most here
  p <- sort(pmax(p, .Machine$double.xmin))
  n <- length(p)
  idx <- thin_qq_idx(n)
  data.table(
    instrument_cat = group_label,
    expected = -log10((idx - 0.5) / n),
    observed = -log10(p[idx])
  )
}

df_qq <- rbindlist(list(
  make_qq_df(df_msmr_tenk10k[instrument_cat == "Single instrument", p_SMR_multi], "Single instrument"),
  make_qq_df(df_msmr_tenk10k[instrument_cat == "Multiple instruments (>1)", p_SMR_multi], "Multiple instruments (>1)")
))
df_qq[, instrument_cat := factor(instrument_cat, levels = levels(df_msmr_tenk10k$instrument_cat))]

(p_qq <- ggplot(df_qq, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = instrument_cat), size = 1.1, alpha = 0.6) +
  scale_color_manual(
    values = instrument_color, name = NULL,
    guide = guide_legend(override.aes = list(shape = "square", size = 3, alpha = 1))
  ) +
  labs(
    x = expression(Expected ~ -log[10](italic(P)[MR])),
    y = expression(Observed ~ -log[10](italic(P)[MR])),
    # title = "QQ plot of the MR test p-value by number of instruments",
    tag = "c"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.01,0.99),
        panel.grid.major.y = element_line(color = "gray95"),
        legend.justification = c(0, 1))
)

# ---------------------------------------------------------------------------
# Panel d: sensitivity-test detection rate (HEIDI, MRLink2 pleiotropy) among
#    MR-significant hits (mr == TRUE), by number of instruments - reuses
#    panel a's instrument_n_bin so the two panels share an x-axis. Cochran's Q
#    is not shown here since it is structurally unavailable for
#    single-instrument tests (see panel b caveat, above), which would leave
#    its detection rate undefined at n = 1; HEIDI and MRLink2 pleiotropy are
#    both computable regardless of instrument count, so a rising detection
#    rate here is genuine evidence that pleiotropy is easier to detect with
#    more instruments, not an artifact of test availability.
# ---------------------------------------------------------------------------
df_mr_hits <- df_msmr_tenk10k[mr == TRUE]
df_mr_hits[, `:=`(
  heidi_rej = p_HEIDI < 0.05,
  mrlink2_rej = psigmay_mrlink2 < 0.05,
  heidi_test = !is.na(p_HEIDI),
  mrlink2_test = !is.na(psigmay_mrlink2)
)]

df_sensitivity_rate <- df_mr_hits[, .(
  `HEIDI` = sum(heidi_rej, na.rm = TRUE) / sum(heidi_test, na.rm = TRUE),
  `MR-link-2` = sum(mrlink2_rej, na.rm = TRUE) / sum(mrlink2_test, na.rm = TRUE)
), by = instrument_n_bin] |>
  pivot_longer(-instrument_n_bin, names_to = "test", values_to = "rejection_rate")

test_color <- paletteer_d("ggthemes::Tableau_10")[3:4] |>
  set_names(c("HEIDI", "MR-link-2"))

(p_panel_d <- ggplot(df_sensitivity_rate, aes(x = instrument_n_bin, y = rejection_rate, color = test, group = test)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = test_color, name = NULL) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Number of instruments",
    y = "Rejection rate (P < 0.05)",
    title = str_wrap("Pleiotropy detection among MR-significant hits", width = 28),
    tag = "d"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "inside",
        legend.justification = c(1, 0),
        legend.position.inside = c(0.99, 0.05),
        panel.grid.major.y = element_line(color = "gray95"),
        plot.title = element_text(size = 11),
        # half-width panel needs a positive gap here - the shared -10pt pull
        # applied to full-width panels below would collide with this panel's
        # (longer, %-formatted) axis text, so this is set again after
        # p_combined's shared theme (which would otherwise override it)
        axis.title.y = element_text(margin = margin(r = 5, unit = "points")))
)

# applied per-panel via `+` on each individual plot (rather than shared across
# the whole patchwork via `&`) so panel d's own axis.title.y override, set
# above, isn't clobbered - `&` re-applies its theme to every subplot last,
# which would win over panel d's local override
shared_theme <- theme(plot.tag = element_text(face = "bold", size = 11),
                       axis.title.y = element_text(margin = margin(r = -10, unit = "points")),
                       plot.tag.position = c(0, 1))

p_panel_a <- p_panel_a + shared_theme
p_evidence_n <- p_evidence_n + shared_theme
p_evidence_prop <- p_evidence_prop + shared_theme
p_qq <- p_qq + shared_theme
p_panel_d <- p_panel_d + shared_theme +
  theme(axis.title.y = element_text(margin = margin(r = 5, unit = "points")))

p_combined <- p_panel_a / (p_evidence_n + p_evidence_prop) / (p_qq + p_panel_d) +
  plot_layout(heights = c(1, 1, 1), guides = "keep")

# Extended Data Fig. 3a-d
ggsave(
  "figures/fig-e3.png",
  p_combined,
  device = ragg::agg_png, width = 8, height = 10, units = "in", res = 300,
  scaling = 1.2
)
# cairo_pdf has no `scaling` argument; scale = 1 / 1.2 keeps the same text-to-canvas ratio
ggsave(
  "figures/fig-e3.pdf",
  p_combined,
  device = cairo_pdf, width = 8, height = 10, units = "in",
  scale = 1 / 1.2
)

print(df_n_tested)
print(df_instrument_dist)
print(df_evidence_n)
print(df_sensitivity_rate)
