# Extended Data Figure 9: GPX1 at the Crohn's disease locus
#   a) locus zoom of the Crohn's GWAS and the CD4 TCM and cDC2 eQTLs
#      (locus object from 1-generate_locus_zoom_object.R)
#   b, c) eQTL vs GWAS effects of the MR instruments in CD4 TCM and cDC2
#      (cross plots, from tenk10k-causal scripts/crohns_publication/cross_plot.R)
# Assembled as in tenk10k-causal 2-plot_locus_zoom_fig.R: panel a is locuszoomr's
# base-graphics plot, recorded and placed beside the ggplot cross plots by patchwork
# (which needs gridGraphics for base graphics). Saved at the published size,
# 180 x 195 mm.
library(tidyverse)
library(data.table)
library(patchwork)
library(locuszoomr)
library(glue)
library(ragg)

target_gene = "GPX1"
probe <- "ENSG00000233276"
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

# --- a) locus zoom (as in tenk10k-causal 2-plot_locus_zoom_fig.R) -------------

max_y = 25
draw_locus <- function() {
  multi_layout(nrow = 3,
               plots = {
                 locus_plot(loclist[["GWAS"]], use_layout = FALSE, legend_pos = 'topleft', pcutoff = 5e-08, labels = tenk10k_all, highlight = target_gene, main = final_plot_names[1], ylim = c(0, max_y),
                            cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.6, filter_gene_biotype = "protein_coding", label_x = c(4, -5))
                 locus_plot(loclist[["CD4_TCM"]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = tenk_ct1, highlight = target_gene, main = final_plot_names[2], ylim = c(0, max_y),
                            cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.6, filter_gene_biotype = "protein_coding")
                 locus_plot(loclist[["cDC2"]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = tenk_ct2, highlight = target_gene, main = final_plot_names[3], ylim = c(0, max_y),
                            cex.main = 0.7, cex = 0.6, cex.axis = 0.7, cex.lab = 0.7, cex.text = 0.6, filter_gene_biotype = "protein_coding")
               })
}

# record the plot (as upstream recordPlot()); under Rscript the display list of a
# new device is off, so enable it before drawing
pdf(NULL, width = 90 / 25.4, height = 200 / 25.4)
dev.control("enable")
draw_locus()
c <- recordPlot()
dev.off()

# --- b, c) cross plots ----------------------------------------------------------

df_ins <- arrow::read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.snps4msmr.parquet.gz") |>
  filter(phenotype == "crohns")

make_plot_mr <- function(
    df_x, df_y, x_ins,
    se_mult = 1,
    label = list(
      x = bquote(beta[eQTL]),
      y = bquote(beta[GWAS]),
      lead_var = "Lead instrument",
      clump_var = "MR instrument"
    ),
    col_scheme = list(
      lead_var = "#6A1B9A",
      clump_var = "red3",
      other_var = "gray75"
    )) {

  df_plot <- inner_join(df_x, df_y, by = c("variant_id", "chr", "pos")) %>%
    # align beta for ea.x
    mutate(ins = variant_id %in% x_ins,
           b_eax.x = b.x,
           b_eax.y = ifelse(ea.x == ea.y, b.y, -b.y),
           li_eax.x = b_eax.x - se_mult*se.x,
           ui_eax.x = b_eax.x + se_mult*se.x,
           li_eax.y = b_eax.y - se_mult*se.y,
           ui_eax.y = b_eax.y + se_mult*se.y,
           lead_var = p.x == min(p.x)
    )

  df_lead_var <- df_plot[ins == TRUE & lead_var == TRUE]
  b_smr <- df_lead_var[,b_eax.y/b_eax.x]
  xlim <- c(min(c(df_plot[ins == TRUE, li_eax.x], df_plot[,b_eax.x])),
            max(c(df_plot[ins == TRUE, ui_eax.x], df_plot[,b_eax.x])))
  ylim <- c(min(c(df_plot[ins == TRUE, li_eax.y], df_plot[,b_eax.y])),
            max(c(df_plot[ins == TRUE, ui_eax.y], df_plot[,b_eax.y])))

  b_smr_label <- list(bquote(beta["MR lead"] == .(scales::number(b_smr, 0.01))))

  ggplot(df_plot, aes(x = b_eax.x, y = b_eax.y, ymin = li_eax.y, ymax = ui_eax.y,
                           xmin = li_eax.x, xmax = ui_eax.x)) +
    theme_classic(base_size = 7) +
    geom_point(color = "gray", alpha = 0.5, size = 0.6) +
    geom_abline(slope = b_smr, intercept = 0, color = col_scheme$lead_var, linewidth = 0.4,
               linetype = "solid") +
    geom_point(aes(color = col_scheme$lead_var), data = df_lead_var, fill = "white",
               shape = "circle filled", size = 2.5, stroke = 0.7) +
    geom_errorbar(data = ~filter(.x, ins == TRUE), width = 0.002, linewidth = 0.3, color = col_scheme$clump_var) +
    geom_errorbar(data = ~filter(.x, ins == TRUE), orientation = "y", width = 0.01, linewidth = 0.3, color = col_scheme$clump_var) +
    geom_point(aes(color = col_scheme$clump_var), data = ~filter(.x, ins == TRUE),
               shape = "circle", size = 1.2, stroke = 0.4) +
    scale_color_identity(labels = c(label$clump_var, label$lead_var),
                         breaks = c(col_scheme$clump_var, col_scheme$lead_var),
                         name = NULL,
                         guide = guide_legend()) +
    geom_hline(yintercept = 0, color = "gray20", linewidth = 0.3, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "gray20", linewidth = 0.3, linetype = "dashed") +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE, clip = "off") +
    annotate("label", label = b_smr_label, parse = TRUE,
             # x = xlim[2]*1.1 , y = xlim[2]*1.1 * b_smr,
             x = if (b_smr > 0) xlim[1] - 0.01 *xlim[1] else xlim[2] - 0.01 * xlim[2],
             y = ylim[2] - 0.01 * ylim[2],
             hjust = if (b_smr > 0) 0 else 1,
             vjust = 1,
             size = 7/.pt, color = col_scheme$lead_var, fill = "white",
             linewidth = 0,
             fontface = "bold", lineheight = 1.2) +
    labs(x = label$x, y = label$y) +
    theme(axis.ticks.length = unit(0.25, "lines"),
          plot.title = element_text(size = 7),
          legend.position = "inside",
          legend.key.size = unit(0.8, "lines"),
          legend.box.background = element_rect(fill = NA, color = NA),
          legend.position.inside = if(b_smr < 0) c(0.01,0.01) else c(0.99,0.01),
          legend.justification = if(b_smr < 0) c(0,0) else c(1,0),
          panel.grid.minor = element_blank())
}

df_metadata <- tibble(
  pheno = c("crohns", "crohns"),
  biosample = c("CD4_TCM", "cDC2"),
  probe = c(probe, probe),
  gene = c(target_gene, target_gene)
) |>
  mutate(
    file_x = glue("results/smr_locus/tenk10k_phase1/{biosample}/{pheno}/{probe}.eqtl.tsv"),
    file_y = glue("results/smr_locus/tenk10k_phase1/{biosample}/{pheno}/{probe}.gwas.tsv"),
    df_x = map(file_x, fread),
    df_y = map(file_y, fread),
    x_ins = pmap(list(pheno, biosample, probe), function(x, y, z) {
      df_ins |>
        filter(phenotype == x, biosample == y, gene == z) |>
        pull(variant)
    }),
    plot = pmap(list(df_x, df_y, x_ins, gene, biosample), function(a, b, c, d, e) {
      make_plot_mr(a,b , c) +
        labs(title = bquote(italic(.(d)) ~ "in" ~.(str_replace_all(e, "_", " ")))) +
        theme(plot.margin = margin(5.5, 10, 5.5, 5.5))  # room for the last x-axis label
    })
  )

# --- assemble (as upstream) ------------------------------------------------------

a <- df_metadata$plot[[1]]
b <- df_metadata$plot[[2]]

combined <- wrap_elements(panel = ~replayPlot(c)) + (a / b) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold"))) +
  plot_layout(widths = c(1.3, 1))

# Extended Data Figure 9
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/fig-e9.pdf", combined, width = 180, height = 195, units = "mm", device = cairo_pdf, bg = "white")
ggsave("figures/fig-e9.png", combined, width = 180, height = 195, units = "mm", dpi = 500, device = agg_png, bg = "white")
