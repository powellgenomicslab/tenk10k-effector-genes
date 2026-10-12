# locus zoom crohns
# library(rtracklayer)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
library(ragg)
library(fs)
library(scales)

# cross plot MR
# interactive test

make_plot_locus <- function(
    df_assoc, df_ld, x_ins, lead_var, lead_var_pos,
    chr, start, end,
    label = list(
      x = paste0("Chromosome ", chr, " position (Mbp)"),
      y = bquote(-log[10]~italic(P)),
      lead_var = "Lead instrument",
      clump_var = "MR instrument"
    ),
    col_scheme = list(
      lead_var = "#6A1B9A",
      clump_var = "red3",
      other_var = "gray75",
      rsq_cat = rev(c("#E65100", "#F9A825", "#76FF03", "#18FFFF", "#1A237E")),
      rsq_missing = "gray75"
    )) {
  
  df_assoc[df_ld, r2_top_var := i.R2, on = c("variant_id" = "SNP_B")]
  df_assoc[, ins := variant_id %in% x_ins]
  
  df_assoc[order(r2_top_var, na.last = FALSE)] %>% 
      ggplot(aes(x = pos, y = -log10(p))) +
      theme_minimal() +
      # plot non-lead SNPs
      geom_point(aes(color = r2_top_var, 
                     fill = after_scale(alpha(color, 0.6))),
                 data = ~.x[ins == FALSE],
                 shape = "circle filled",
                 size = 1.5, stroke = 0.5) +
      # plot lead SNPs
      geom_vline(xintercept = lead_var_pos, color = col_scheme$lead_var, linewidth = 0.5,
                 linetype = "dashed") + 
      geom_point(aes(fill = "white"),
                 shape = "circle filled",
                 data = ~.x[variant_id == lead_var], stroke = 1,
                 color = col_scheme$lead_var, size = 3.5) +
      geom_point(aes(fill = col_scheme$clump_var),
                 shape = "circle filled",
                 data = ~.x[ins == TRUE], alpha = 1,
                 size = 3, stroke = 0) +
    scale_fill_identity(labels = c(label$clump_var, label$lead_var),
                        breaks = c(col_scheme$clump_var, "white"),
                         name = NULL) +
    scale_color_stepsn(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        right = FALSE,
        colours = col_scheme$rsq_cat,
        na.value = col_scheme$rsq_missing,
        name = bquote("LD"~italic(R) ^ 2),
        guide = guide_coloursteps(order = 2, frame.linewidth = 0.5, ticks.linewidth = 0.5,
                                  frame.colour = "black", ticks = TRUE, ticks.colour = "black")
      ) +
      scale_x_continuous(labels = label_number(0.01, scale = 1e-6),
                         limits = c(start, end),
                         expand = expansion(mult = 0.01)) +
      labs(y = label$y, x = label$x) +
      coord_cartesian(clip = "off") +
      theme(legend.spacing = unit(0, "lines"),
            legend.key.width = unit(1.1, "lines"),
            legend.key.height = unit(0.8, "lines"),
            legend.title = element_text(size = 10, face = "bold"),
            axis.title.x = element_blank(),
            panel.grid.minor = element_blank())
}

# SMREffectPlot(data_smr)
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
  
  df_lead_var <- df_plot[p.x == min(p.x)]
  b_smr <- df_lead_var[,b_eax.y/b_eax.x]
  xlim <- c(min(c(df_plot[ins == TRUE, li_eax.x], df_plot[,b_eax.x])),
            max(c(df_plot[ins == TRUE, ui_eax.x], df_plot[,b_eax.x])))
  ylim <- c(min(c(df_plot[ins == TRUE, li_eax.y], df_plot[,b_eax.y])),
            max(c(df_plot[ins == TRUE, ui_eax.y], df_plot[,b_eax.y])))
            
  b_smr_label <- list(bquote(beta["MR"] == .(scales::number(b_smr, 0.01))))
  
  ggplot(df_plot, aes(x = b_eax.x, y = b_eax.y, ymin = li_eax.y, ymax = ui_eax.y,
                           xmin = li_eax.x, xmax = ui_eax.x)) +
    theme_classic() +
    geom_point(color = "gray", alpha = 0.5) +
    geom_abline(slope = b_smr, intercept = 0, color = col_scheme$lead_var, linewidth = 0.5,
               linetype = "solid") +
    # geom_errorbar(data = df_lead_var, width = 0.01, color = col_scheme$lead_var) +
    # geom_errorbarh(data = df_lead_var, width = 0.01, color = col_scheme$lead_var) +
    geom_point(aes(color = col_scheme$lead_var), data = df_lead_var, fill = "white",
               shape = "circle filled", size = 4, stroke = 1) +
    geom_errorbar(data = ~filter(.x, ins == TRUE), width = 0.01, color = col_scheme$clump_var) +
    geom_errorbarh(data = ~filter(.x, ins == TRUE), width = 0.01, color = col_scheme$clump_var) +
    geom_point(aes(color = col_scheme$clump_var), data = ~filter(.x, ins == TRUE),
               shape = "circle", size = 2, stroke = 0.5) +
    scale_color_identity(labels = c(label$clump_var, label$lead_var),
                         breaks = c(col_scheme$clump_var, col_scheme$lead_var),
                         name = NULL,
                         guide = guide_legend()) +
    geom_hline(yintercept = 0, color = "gray20", linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "gray20", linewidth = 0.5, linetype = "dashed") +
    coord_fixed(xlim = xlim, ylim = ylim, expand = TRUE, clip = "off") +
    annotate("label", label = b_smr_label, parse = TRUE,
             # x = xlim[2]*1.1 , y = xlim[2]*1.1 * b_smr,
             x = if (b_smr > 0) xlim[1] - 0.01 *xlim[1] else xlim[2] - 0.01 * xlim[2],
             y = ylim[2] - 0.01 * ylim[2],
             hjust = if (b_smr > 0) 0 else 1,
             vjust = 1,
             size = 9/.pt, color = col_scheme$lead_var, fill = "white",
             linewidth = 0,
             fontface = "bold", lineheight = 1.2) +
    labs(x = label$x, y = label$y) +
    theme(axis.ticks.length = unit(0.25, "lines"),
          legend.position = "inside",
          legend.box.background = element_rect(fill = NA, color = NA),
          legend.position.inside = if(b_smr < 0) c(0.01,0) else c(0.01,1),
          legend.justification = if(b_smr < 0) c(0,0) else c(0,1),
          # plot.margin = margin(r = 5, unit = "lines"),
          panel.grid.minor = element_blank())
}


combine_locus_plot <- function(pheno, biosample, probe,
                               label = list(
                                 gwas = NULL,
                                 cell = NULL,
                                 gene = NULL
                               ),
                               start_win_bp = 1e5, end_win_bp = 1e5,
                               base_dir = "results/smr_locus/tenk10k_phase1") {
  df_gwas <- fread(fs::path(base_dir, biosample, pheno, paste0(probe, ".gwas.tsv")))
  df_eqtl <- fread(fs::path(base_dir, biosample, pheno, paste0(probe, ".eqtl.tsv")))
  df_ld <- fread(fs::path(base_dir, biosample, pheno, paste0(probe, ".ld")))
  eqtl_ins <- fread(fs::path(base_dir, biosample, pheno, paste0(probe, ".clump")))$SNP
  
  vars <- intersect(df_gwas$variant_id, df_eqtl$variant_id)
  df_gwas <- df_gwas[variant_id %in% vars]
  df_eqtl <- df_eqtl[variant_id %in% vars]
  
  p_thresh <- as.numeric(readLines(fs::path("resources/smr/tenk10k_phase1", biosample, "pthresh.txt")))
  lead_eqtl <- df_eqtl[include == TRUE & p == min(p), variant_id]
  lead_eqtl_pos <- df_eqtl[variant_id == lead_eqtl,pos]
  chr <- df_gene_annot[ensembl_gene_id == probe, chr]
  start <- df_gene_annot[ensembl_gene_id == probe, start - 1e5]
  end <- df_gene_annot[ensembl_gene_id == probe, end + 1e5]
  
  p_gwas <- make_plot_locus(df_gwas, df_ld, eqtl_ins, lead_eqtl, lead_eqtl_pos, chr, start, end) +
    annotate("label", label = label$gwas, linewidth = 0, fill = alpha("white", 0.5),
             hjust = 0, vjust = 1, x = -Inf, y = Inf, fontface = "bold",
             size = 8/.pt) +
    labs(y = bquote(-log[10]~italic(P)["GWAS"])) +
    theme(axis.text.x = element_blank())
  
  p_eqtl <- make_plot_locus(df_eqtl, df_ld, eqtl_ins, lead_eqtl, lead_eqtl_pos, chr, start, end) + 
    annotate("label", label = paste("\nin", label$cell), linewidth = 0, fill = alpha("white", 0.5),
             hjust = 0, vjust = 1, x = -Inf, y = Inf, fontface = "bold", lineheight = 1.25,
             size = 7/.pt) +
    annotate("label", label = list(bquote(bolditalic(.(label$gene)))),
             linewidth = 0, fill = NA,
             hjust = 0, vjust = 1, x = -Inf, y = Inf, parse = TRUE,
             size = 8/.pt) +
    geom_hline(
      yintercept = -log10(p_thresh), linetype = "dashed", color = "red3",
      linewidth = 0.5, alpha = 1
    ) +
    labs(y = bquote(-log[10]~italic(P)["eQTL"])) +
    annotate("label", label = list(bquote(italic(P)["eQTL"]~"threshold")), 
             parse = TRUE, x = -Inf, y = -log10(p_thresh), hjust = 0,
             linewidth = 0,
             color = "red3", size = 3, fill = alpha("white", 0.5), vjust = -0.1) +
    theme(axis.ticks.x = element_line(),
          axis.ticks.length.x = unit(0.25, "lines"),
          axis.title.x = element_text(size = 9))
  
  p_mr <- make_plot_mr(df_eqtl, df_gwas, eqtl_ins) +
    theme(plot.margin = margin(r = 0.5, b = 0.4, unit = "lines"),
          aspect.ratio = 1,
          legend.position = "bottom")
  
  locus_plot <- ((p_gwas / p_eqtl) &
                   guides(fill = "none") &
                   theme(plot.margin = margin(),
                         legend.key.size = unit(0.85, "lines"),
                         panel.border = element_rect(linewidth = 0.3),
                         axis.ticks.y = element_line(),
                         axis.title = element_text(size = 9),
                         axis.ticks.length.y = unit(0.25, "lines"))) +
    plot_layout(ncol = 1)
  
  (p_mr + locus_plot) +
    plot_layout(widths = c(0.45, 0.55))
}

df_targets <- readxl::read_xlsx("resources/misc/crohns_example_gene.xlsx") %>% 
  filter(!is.na(order)) %>% 
  arrange(order) %>% 
  group_by(pheno, biosample, probe, order) %>% 
  nest() %>% 
  mutate(label = map(data, as.list)) %>%
  mutate(
    plots = pmap(list(pheno, biosample, probe, label), combine_locus_plot)
  )

(mr_example_plots <- wrap_plots(df_targets$plots, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom",
        legend.title.position = "left",
        legend.title = element_text(size = 9, vjust = 1),
        legend.box.spacing = unit(0, "lines"),
        legend.box.margin = margin(),
        legend.spacing.x = unit(5, "lines"),
        legend.margin = margin(t = 0.5, unit = "lines"),
        axis.title = element_text(size = 9)
        # plot.margin = margin()
        ))
