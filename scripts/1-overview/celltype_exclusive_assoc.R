# Supplementary Figure 4: cell-type-exclusive MR associations per cell type
# Extracted unchanged from scripts/publication/mr_example/4_pub_combined.R (tenk10k-causal),
# where it was appended to the Fig 3 script.

library(tidyverse)
library(data.table)
library(paletteer)
library(ragg)

source("scripts/0-preprocess/preprocess_results.R")

# unique cell type to address reviewer's comment
df_msmr <- df_msmr_tenk10k[mr == TRUE]
group_col <- expr(magma_gene)
df_tally_unique <- df_msmr %>% 
    group_by(probeID, phenotype) %>%
    filter(n() == 1) |> 
    group_by(cell_type, pheno_cat) %>%
    tally(name = "n") %>% 
    arrange(n)

# bar chart of N by cell type
p_celltype_specific <- ggplot(df_tally_unique, aes(x = fct_reorder(cell_type, n, sum), y = n)) +
    geom_col(aes(fill = pheno_cat)) +
    theme_classic() +
    labs(y = "N cell-type specific MR associations",
         x = NULL) +
    # annotate with total number
    geom_text(aes(label = n),
              data = df_tally_unique %>% group_by(cell_type) %>% summarise(n = sum(n)),
              nudge_y = 50, angle = 90, size = 9/.pt, hjust = 0) +
    scale_fill_paletteer_d("ggthemes::Tableau_10", name = "Phenotype category") +
    coord_cartesian(clip = "off") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(1,1,1,1, unit = "lines"))

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig-s4.png", p_celltype_specific, width = 7, height = 4,
       device = ragg::agg_png, scaling = 0.8)
ggsave("figures/fig-s4.pdf", p_celltype_specific, width = 7 / 0.8, height = 4 / 0.8,
       device = cairo_pdf)
