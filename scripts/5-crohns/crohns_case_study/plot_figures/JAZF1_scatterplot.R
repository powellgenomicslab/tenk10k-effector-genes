library(tidyverse)
library(locuszoomr)
library(ggrepel)
# scatterplot firections of effect 

plot_data <- readRDS("resources/crohns_case_study/figures/JAZF1_data.RDS")

join <- full_join(plot_data[["CD4_Naive"]]$data, plot_data[["cDC2"]]$data, by = "MarkerID", suffix = c("_CD4Naive", "_cDC2"))

# Get actual intrsuments from resources/sensitivity/celltype/phenotype/look up gene
cd4_instruments <- c("7:28135367:T:C", "7:28009514:C:T", "7:28058805:A:G", "7:28075634:T:C", "7:28122072:G:A", "7:28135367:T:C", "7:28139777:C:T", "7:28185640:G:A", "7:28221043:T:C", "7:28227901:G:C")
cdc2_instruments <-	c("7:28135367:T:C")
plot_data <- join %>% filter(p_CD4Naive < 0.05 | p_cDC2 < 0.05) %>% 
  mutate(marker_label = ifelse(MarkerID %in% c(cd4_instruments, cdc2_instruments), MarkerID, ""))
  
(p <- ggplot(plot_data, aes(x =p_transform_CD4Naive, y = p_transform_cDC2)) + 
    geom_point(aes(fill = ld_CD4Naive), shape = "circle filled", size = 2) +
    geom_label_repel(
      aes(label = marker_label),
      min.segment.length = 0,
      max.overlaps = Inf,
      segment.size = 0.1) +
    geom_point(data = subset(plot_data, MarkerID == "7:28135367:T:C"), 
            shape = "diamond filled", fill = "red3", size = 5) +
  scale_fill_viridis_b(n.breaks = 7) +
  theme_bw() +
  geom_hline(yintercept = 0 , linetype = "dashed") +
  geom_vline(xintercept = 0 , linetype = "dashed") +
  labs(x = "CD4 Naive signed log10 P eQTL", y = "cDC2 signed log10 P eQTL", fill = "LD r2")
)

ggsave("resources/crohns_case_study/figures/JAZF1_signed_eqtl_scatter.png", p, device = ragg::agg_png(),
       width = 5.0, height = 4.5, bg = "white", scaling = 0.7, dpi = 300)

################### 



png(paste0(figs, gene, "_compare_discordant_bsmr_2.png"), width = 5, height = 8, units = "in", res = 300, type = "cairo") # Note it's res, not dpi. 
oldpar <- set_layers(1)
scatter_plot(loclist[[1]], xticks = FALSE, pcutoff = 5e-08, ylim = c(-100, 100), scheme = "blue", showLD = FALSE)
scatter_plot(loclist[[2]], xticks = FALSE, pcutoff = 5e-08, scheme = "orange", showLD = FALSE, pch = 22, add = TRUE)
genetracks(loclist[[2]])
par(oldpar)
dev.off()





