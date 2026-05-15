library(ggplot2)
library(tidyverse)
library(patchwork)

mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS")
mr_all <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_all.RDS")
cell_map <- read_tsv("resources/metadata/cell_map.tsv")
colours <- cell_map %>% select(cell_type, color) %>% deframe()


# sig
gene_spec <- table(mr$cell_type) %>% 
  as.data.frame() %>% 
  rename(cell_type = Var1, Gene_Number = Freq) %>% 
  arrange(Gene_Number)

# all 
gene_spec_all <- table(mr_all$cell_type) %>% 
  as.data.frame() %>% 
  rename(cell_type = Var1, Gene_Number = Freq) %>% 
  arrange(Gene_Number)


# create a bar plot of the number of sig/causal genes per cell type 

(histogram_of_gene_numbers <- ggplot(gene_spec, aes(x = fct_reorder(cell_type, Gene_Number), 
                                                    y = Gene_Number)) + 
    geom_col(aes(fill = cell_type), color = "#e9ecef", alpha=0.9) +
    # geom_col(fill = "#69b3a2", colour = "black", alpha=0.9) +
    ggtitle("Crohn's Disease TenK10K single-cell MR significant genes per cell type") +
    # labs(x = "Cell Types", y = "Number of causal eGenes") +
    scale_y_continuous(n.breaks = 8) +
    # scale_x_continuous(breaks = 1:max(gene_spec$Gene_Number)) +
    scale_fill_manual(values = colours) +
    geom_text(aes(label = Gene_Number), vjust = -0.5, size = 2) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=12),
      panel.grid.major.y = element_line(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
)


# create a bar plot of the number of tested genes per cell type 

(histogram_of_gene_numbers_all <- ggplot(gene_spec_all, aes(x = fct_reorder(cell_type, Gene_Number), 
                                                            y = Gene_Number)) + 
    geom_col(aes(fill = cell_type), color = "#e9ecef", alpha=0.9) +
    # geom_col(fill = "#69b3a2", color = "black", alpha=0.9) +
    ggtitle("Number of genes tested per cell type in TenK10K single-cell MR") +
    # labs(x = "Cell Types", y = "Number of eGenes") +
    scale_y_continuous(n.breaks = 7) +
    # scale_x_continuous(breaks = 1:max(gene_spec$Gene_Number)) +
    scale_fill_manual(values = colours) +
    geom_text(aes(label = Gene_Number), vjust = -0.5, size = 2) +
    coord_cartesian(clip = "off") +
    # coord_flip()+
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=12),
      panel.grid.major.y = element_line(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
)

(p <- histogram_of_gene_numbers + histogram_of_gene_numbers_all + plot_layout(nrow = 2))

FIG3 = (histogram_of_gene_numbers / histogram_of_gene_numbers_all) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))

ggsave("resources/crohns_case_study/figures/gene_numbers_across_celltypes_corrected.png", FIG3, device = ragg::agg_png,
       width = 12, height = 8, bg = "white", scaling = 1.2, dpi = 300)

