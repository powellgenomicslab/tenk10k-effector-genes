library(ggplot2)
library(tidyverse)
library(patchwork)

mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS")
eqtlgen <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/eqtlgen_crohns_sig_genes.RDS")
cell_map <- read_tsv("resources/metadata/cell_map.tsv")
colours <- cell_map %>% select(cell_type, color) %>% deframe()

ct_spec <- table(mr$probeID) %>% 
  as.data.frame() %>% 
  rename(probeID = Var1, ct_number = Freq) %>% 
  arrange(ct_number) %>% 
  mutate(cell_type_specific = ifelse(ct_number == 1, T, F)) %>% 
  mutate(eqtlgen_sig = ifelse(probeID %in% eqtlgen$probeID , "TenK10K + eQTLGen", "TenK10K only")) %>% 
  mutate(magma_sig = ifelse(probeID %in% mr[mr$magma_gene == T,]$probeID , "MR + GWAS", "MR only")) %>% 
  group_by(ct_number) %>% 
  mutate(percentage = (sum(magma_sig == "MR only"))/n()*100) %>% 
  ungroup()
  
  
# this plot includes the eqtl gen significance
(gwas_plot <- ggplot(ct_spec, aes(x = ct_number, fill = magma_sig)) + 
  # geom_bar(fill = "#69b3a2", color = "#e9ecef", alpha=0.9) +
  geom_bar(colour = "black", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("#4C6C94", "grey")) +
    
  # scale_fill_brewer(palette = "Pastel1") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "N cell types with MR associations", y = "N gene-trait MR associations", fill = "MAGMA", title = "TenK10K single-cell MR vs. GWAS for Crohn's Disease") +
  scale_y_continuous(n.breaks = 7) +
  scale_x_continuous(breaks = 1:max(ct_spec$ct_number), n.breaks = 2) +
  #geom_line(aes(x=ct_number, y = percentage)) +
  # coord_flip() +
  # geom_bar(filter(eqtlgen_sig == T), fill = "red")) + 
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8,0.7),
    legend.title = element_blank(),
    plot.title = element_text(size=14),
    panel.grid.major.y = element_line(),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

)

(eqtlgen_plot <- ggplot(ct_spec, aes(x = ct_number, fill = eqtlgen_sig)) + 
    # geom_bar(fill = "#69b3a2", color = "#e9ecef", alpha=0.9) +
    geom_bar(colour = "black", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c("#4C6C94", "grey")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = "N cell types with MR associations", y = "N gene-trait MR associations", fill = "eqtlgen", title =  "TenK10K single-cell MR vs. eQTLgen bulk whole blood for Crohn's Disease") +
    scale_y_continuous(n.breaks = 7) +
    scale_x_continuous(breaks = 1:max(ct_spec$ct_number)) +
    # coord_flip() +
    # geom_bar(filter(eqtlgen_sig == T), fill = "red")) + 
    theme_classic() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.8,0.7),
      legend.title = element_blank(),
      plot.title = element_text(size=14),
      panel.grid.major.y = element_line(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    )
  
)

FIG3 = (gwas_plot / eqtlgen_plot) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))


ggsave("resources/crohns_case_study/figures/magma_bulkMR_comparison_by_celltype_numbers.png", FIG3, device = ragg::agg_png,
       width = 10.0, height = 12, bg = "white", scaling = 1.2, dpi = 300)

