library(tidyverse)
library(cowplot)
library(locuszoomr)

target_gene = "GPX1"
loclist <- readRDS(paste0("resources/crohns_case_study/locus_zoom_data/", target_gene, "/tenk10k_", target_gene, "_locus_plot_data.RDS"))

# Clean the position columns in your locus objects before plotting
for (name in names(loclist)) {
  loclist[[name]]$data$pos <- as.numeric(as.character(loclist[[name]]$data$pos))
}

final_plot_names = c("Crohns Disease GWAS", 
                     "CD4 TCM TenK10K eQTL",
                     "cDC2 TenK10K eQTL")

# Get instruments from resources/sensitivity/smr/{CELLTYPE}/{trait}/allchr.snps4msmr.list
tenk_ct1 = "3:49369275:G:A"
tenk_ct2 = "3:49361826:G:A"
tenk10k_all <- unique(c(tenk_ct1, tenk_ct2))

# Add these parameters to your locus_plot() calls:
#   • cex = 0.8 (Scales down the size of the scatter plot data points)
# • cex.axis = 0.7 (Shrinks the tick numbers on the X and Y axes)
# • cex.lab = 0.8 (Shrinks axis titles like Chromosome position or -log10(p-value))
# • cex.text = 0.6 (Shrinks the gene name annotations at the bottom)
# label_x = c(4,  -5),  # 1st label moves 4% RIGHT,  2nd label moves 5% LEFT
# label_y = c(2,   5)   # 1st label moves 2% UP,     2nd label moves 5% UP

# Extended Data Figure 9
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
png("figures/fig-e9.png", width = 90, height = 170, units = "mm", res = 500)
max_y = 25
# max_y_2 = 70
multi_layout(nrow = 3,
             plots = {
               locus_plot(loclist[["GWAS"]], use_layout = FALSE, legend_pos = 'topleft', pcutoff = 5e-08, labels = tenk10k_all, highlight = target_gene, main = final_plot_names[1], ylim = c(0, max_y),
                          cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.5, filter_gene_biotype = "protein_coding", label_x = c(4, -5))
               locus_plot(loclist[["CD4_TCM"]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = tenk_ct1, highlight = target_gene, main = final_plot_names[2], ylim = c(0, max_y),
                          cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.5, filter_gene_biotype = "protein_coding")
               locus_plot(loclist[["cDC2"]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = tenk_ct2, highlight = target_gene, main = final_plot_names[3], ylim = c(0, max_y),
                          cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.5, filter_gene_biotype = "protein_coding")
                 })
dev.off()

# convert to mm from default inches
# width = 90 / 25.4,
pdf("figures/fig-e9.pdf", width = (90 / 25.4), height = (170 / 25.4))

max_y = 25
# max_y_2 = 70
multi_layout(nrow = 3,
             plots = {
               locus_plot(loclist[["GWAS"]], use_layout = FALSE, legend_pos = 'topleft', pcutoff = 5e-08, labels = tenk10k_all, highlight = target_gene, main = final_plot_names[1], ylim = c(0, max_y),
                          cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.5, filter_gene_biotype = "protein_coding", label_x = c(4, -5))
               locus_plot(loclist[["CD4_TCM"]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = tenk_ct1, highlight = target_gene, main = final_plot_names[2], ylim = c(0, max_y),
                          cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.5, filter_gene_biotype = "protein_coding")
               locus_plot(loclist[["cDC2"]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = tenk_ct2, highlight = target_gene, main = final_plot_names[3], ylim = c(0, max_y),
                          cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.5, filter_gene_biotype = "protein_coding")
             })
dev.off()



