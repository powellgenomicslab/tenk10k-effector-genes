# post-hoc power for the MoA-MR concordance binomial test (supplementary table/figure)

library(tidyverse)
library(data.table)
library(ragg)
library(scales)
library(patchwork)

df_concordance <- fread("results/otp/26.03/otp_concordance_by_phase.tsv")

# exact one-sided binomial power at alpha = 0.05: P(X >= critical | p_true), where critical is the smallest X with P(X >= critical | p = 0.5) <= alpha
binom_power <- function(n, p_true, alpha = 0.05) {
  critical <- qbinom(1 - alpha, n, 0.5) + 1
  1 - pbinom(critical - 1, n, p_true)
}

# project each phase's observed rate to larger sample sizes (up to ~ the full T-I pair count, n = 105), holding the rate fixed
multipliers <- c(1, 2, 3, 4, 5)
phase_levels <- c("Phase I", "Phase II", "Phase III", "Approved")

df_power <- df_concordance %>%
  select(phase, n, n_concordant, concordance) %>%
  cross_join(tibble(multiplier = multipliers)) %>%
  mutate(n_proj = round(n * multiplier),
         n_concordant_proj = round(concordance * n_proj),
         concordance_proj = n_concordant_proj / n_proj,
         is_current = multiplier == 1,
         # discrete axis: keeps 1x-5x evenly spaced so labels don't crowd together the way they would on a linear numeric scale
         multiplier_label = factor(paste0(multiplier, "x"), levels = paste0(multipliers, "x")),
         p_binomial = map2_dbl(n_concordant_proj, n_proj,
                                ~binom.test(.x, .y, p = 0.5, alternative = "greater")$p.value),
         power_binomial = binom_power(n_proj, concordance),
         phase = factor(phase, phase_levels))

# fwrite(df_power, "results/otp/26.03/otp_power_projection.tsv", sep = "\t")

# binomial power at increasing projected sample size, current sample size highlighted; x-axis is a multiple of each phase's own observed n so the four panels share one scale
# N concordant / N total at each projected multiple, annotated below the x-axis
df_n_label <- df_power %>%
  distinct(phase, multiplier_label, n_concordant_proj, n_proj) %>%
  mutate(label = paste0(n_concordant_proj, "/", n_proj))

(p_power <- df_power %>%
  mutate(method = "Binomial test", power = power_binomial) %>%
  ggplot(aes(x = multiplier_label, y = power, color = method, group = method)) +
  theme_bw() +
  geom_line() +
  geom_point(aes(shape = is_current, size = is_current)) +
  geom_text(data = df_n_label, aes(x = multiplier_label, y = -Inf, vjust = -0.5, label = label),
            inherit.aes = FALSE, size = 7 / .pt, color = "grey30", family = "Arial Narrow") +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 3.5, `FALSE` = 2), guide = "none") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60") +
  facet_wrap(~phase, nrow = 1) +
  scale_color_manual(values = c("Binomial test" = "#3B7DD8")) +
  scale_y_continuous(labels = percent, breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.2, 0.05))) +
  annotate("text", x = 1, y = -0.05, vjust = 1, hjust = 0,family = "Arial Narrow",
           label = "N concordant / N total", size = 7/.pt, fontface= "bold") +
  labs(x = "Projected sample size",
       y = "Estimated power",
       color = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
)
# ggsave("figures/strict/otp/otp_power_projection.png", p_power,
#        device = agg_png, bg = "white", width = 9, height = 3.8, dpi = 300, scaling = 1.2)

# projected -log10(P) at increasing projected sample size, current sample size highlighted
(p_pvalue <- df_power %>%
  mutate(method = "Binomial test", p_value = p_binomial) %>%
  ggplot(aes(x = multiplier_label, y = -log10(p_value), color = method, group = method)) +
  theme_bw() +
  geom_line() +
  geom_point(aes(shape = is_current, size = is_current)) +
  geom_text(data = df_n_label, aes(x = multiplier_label, y = -Inf, vjust = -0.5, label = label),
            inherit.aes = FALSE, size = 7 / .pt, color = "grey30", family = "Arial Narrow") +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 3.5, `FALSE` = 2), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  facet_wrap(~phase, nrow = 1) +
  scale_color_manual(values = c("Binomial test" = "#3B7DD8")) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05))) +
  annotate("text", x = 1, y = -0.05, vjust = 1, hjust = 0, family = "Arial Narrow",
           label = "N concordant / N total", size = 7/.pt, fontface = "bold") +
  labs(x = "Projected sample size",
       y = expression(-log[10](italic(P))),       
       caption = "▲: observed sample size",
       color = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
)

plots <- (p_power / p_pvalue) +
  plot_layout(guide = "collect") &
  theme(legend.position = "none")

ggsave("figures/revision2/otp/otp_power_projection_pvalue.png", plots,
       device = agg_png, bg = "white", width = 9, height = 6, dpi = 300, scaling = 1.2)

# write to supplementary files

source("scripts/util/write_table.R")
write_table(df_power, "otp_power")
