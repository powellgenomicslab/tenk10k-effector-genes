# post-hoc power for the MoA-MR concordance test (supplementary table/figure)
# holds each phase's observed concordance rate fixed and projects the power
# of the binomial test and the permutation test to detect it at increasing
# sample sizes

library(tidyverse)
library(data.table)
library(ragg)
library(scales)
library(patchwork)

df_concordance <- fread("results/otp/26.03/otp_concordance_by_phase.tsv")

# exact one-sided binomial power at alpha = 0.05: reject H0 (p = 0.5) once
# the concordant count reaches the smallest value with P(X >= critical | p = 0.5) <= alpha;
# power is the probability of reaching that value if the true rate is p_true
binom_power <- function(n, p_true, alpha = 0.05) {
  critical <- qbinom(1 - alpha, n, 0.5) + 1
  1 - pbinom(critical - 1, n, p_true)
}

# Monte Carlo behaviour of the phase-scoped permutation test (no closed
# form): simulate n_sims datasets of size n with a true concordance rate of
# p_true, run the permutation test on each, and summarise the resulting
# p-values as (i) the fraction reaching significance at alpha (power) and
# (ii) the median p-value, for comparison against the exact binomial test
# set.seed(123)
# simulate_permutation_p <- function(n, p_true, n_sims = 500, n_perms = 1000, alpha = 0.05) {
#   p_values <- numeric(n_sims)
#   for (i in seq_len(n_sims)) {
#     observed_dir <- sample(c(1, -1), n, replace = TRUE)
#     is_concordant <- rbinom(n, 1, p_true) == 1
#     predicted_dir <- ifelse(is_concordant, observed_dir, -observed_dir)
#     obs_concordance <- mean(observed_dir == predicted_dir)
#     permuted <- replicate(n_perms, mean(observed_dir == sample(predicted_dir)))
#     p_values[i] <- (sum(permuted >= obs_concordance) + 1) / (n_perms + 1)
#   }
#   list(power = mean(p_values < alpha), median_p = median(p_values))
# }

# project each phase's observed rate to larger sample sizes (up to ~ the
# full T-I pair count, n = 105), holding the rate fixed
multipliers <- c(1, 2, 3, 4, 5)
phase_levels <- c("Phase I", "Phase II", "Phase III", "Approved")

df_power <- df_concordance |>
  select(phase, n, n_concordant, concordance) |>
  cross_join(tibble(multiplier = multipliers)) |>
  mutate(n_proj = round(n * multiplier),
         n_concordant_proj = round(concordance * n_proj),
         concordance_proj = n_concordant_proj / n_proj,
         is_current = multiplier == 1,
         # discrete axis: keeps 1x/2x/5x/10x evenly spaced so labels don't
         # crowd together the way they would on a linear numeric scale
         multiplier_label = factor(paste0(multiplier, "x"), levels = paste0(multipliers, "x")),
         p_binomial = map2_dbl(n_concordant_proj, n_proj,
                                ~binom.test(.x, .y, p = 0.5, alternative = "greater")$p.value),
         power_binomial = binom_power(n_proj, concordance),
         # perm_sim = map2(n_proj, concordance, simulate_permutation_p),
         # power_permutation = map_dbl(perm_sim, "power"),
         # p_permutation = map_dbl(perm_sim, "median_p"),
         phase = factor(phase, phase_levels))
  # select(-perm_sim)

# fwrite(df_power, "results/otp/26.03/otp_power_projection.tsv", sep = "\t")

# supplementary figure: binomial vs. permutation power, current sample size
# highlighted. The x-axis is expressed as a multiple of each phase's own
# observed n (rather than absolute n) so the four panels share one scale.
# N concordant / N total at each projected multiple, annotated below the x-axis
df_n_label <- df_power |>
  distinct(phase, multiplier_label, n_concordant_proj, n_proj) |>
  mutate(label = paste0(n_concordant_proj, "/", n_proj))

(p_power <- df_power |>
  # pivot_longer(c(power_binomial, power_permutation), names_to = "method", values_to = "power") |>
  # mutate(method = recode(method,
  #                         power_binomial = "Binomial test",
  #                         power_permutation = "Permutation test")) |>
  mutate(method = "Binomial test", power = power_binomial) |>
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
  scale_color_manual(values = c("Binomial test" = "#4C6C94FF", "Permutation test" = "#D8553B")) +
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

# supplementary figure: projected -log10(P), current sample size highlighted.
# Permutation P is the median across the simulated datasets above, since a
# single permutation p-value is itself a random draw at a given n and rate.
(p_pvalue <- df_power |>
  # pivot_longer(c(p_binomial, p_permutation), names_to = "method", values_to = "p_value") |>
  # mutate(method = recode(method,
  #                         p_binomial = "Binomial test",
  #                         p_permutation = "Permutation test")) |>
  mutate(method = "Binomial test", p_value = p_binomial) |>
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
  scale_color_manual(values = c("Binomial test" = "#4C6C94FF", "Permutation test" = "#D8553B")) +
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

# Supplementary Figure 16
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig-s16.png", plots,
       device = agg_png, bg = "white", width = 9, height = 6, dpi = 300, scaling = 1.2)
ggsave("figures/fig-s16.pdf", plots,
       device = cairo_pdf, bg = "white", width = 9, height = 6, scale = 1/1.2)

# write to supplementary files

source("scripts/util/write_table.R")
write_table(df_power, "otp_power", 14)  # concordance by phase + power simulation
