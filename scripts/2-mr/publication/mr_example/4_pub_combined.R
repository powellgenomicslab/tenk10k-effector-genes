library(tidyverse)
library(data.table)
library(glue)
library(patchwork)
library(ragg)
library(scales)
library(qvalue)
library(paletteer)
# library(ggnewscale)
library(geomtextpath)
library(arrow)

source("scripts/preprocess_strict.R")

df_msmr <- df_msmr_tenk10k |> 
  filter(mr == TRUE)

gene_universe <- unique(df_msmr_tenk10k$probeID)

df_msmr[df_magma, magma := i.lfdr < 0.05, on = c("probeID" = "GENE", "phenotype" = "phenotype")]
df_msmr[df_msmr_eqtlgen, mr_eqtlgen := i.lfdr_msmr_pheno < 0.05, on = c("probeID", "phenotype")]
df_msmr[is.na(mr_eqtlgen), mr_eqtlgen := FALSE]


mk_plot_tally <- function(group_col, label = list(mr_only = "MR only", mr_other = "MR & MAGMA"),
                          p_title = NULL,
                          col_scheme = list(mr_only = "#A4ABB0FF", mr_other = "#4C6C94FF")){
  df_tally <- df_msmr %>% 
    filter(!is.na({{group_col}})) %>%
    group_by(probeID, phenotype, {{group_col}}) %>%
    tally(name = "n_celltypes") %>%
    group_by(n_celltypes, {{group_col}}) %>%
    tally(name = "n") %>% 
    ungroup() %>% 
    complete(n_celltypes, {{group_col}}, fill = list(n = 0)) %>%  
    mutate(sum_n = sum(n), prop = n / sum_n, .by = n_celltypes) %>% 
    mutate(prop_scaled = rescale(prop, from = c(0, 1), to = range(sum_n)))
  
  ggplot(df_tally, aes(x = n_celltypes, y = n)) +
    geom_col(aes(fill = {{group_col}}), width = 1, linewidth = 0.5,
             position = "stack", color = "black") +
    scale_fill_manual(values = c(col_scheme$mr_only, col_scheme$mr_other),
                      name = NULL,
                      label = c("FALSE" = label$mr_only,
                                "TRUE" = label$mr_other),
                      guide = guide_legend(ncol = 1)) +
    theme_classic() +
    labs(y = "N gene-trait MR associations",
         x = "N cell types with MR associations",
         title = p_title
         ) +
    geom_line(aes(group = {{group_col}}, color = {{group_col}}, y = prop_scaled),
              data = ~filter(.x, !{{group_col}})) +
    geom_textline(aes(group = {{group_col}},
                      label = paste("%", ifelse({{group_col}}, label$mr_other, label$mr_only)),
                      color = {{group_col}}, y = prop_scaled),
                  data = ~filter(.x, !{{group_col}}),
                  size = 9/.pt,
                  vjust = -1, hjust = 0.4) +
    geom_point(aes(color = {{group_col}}, y = prop_scaled), data = ~filter(.x, !{{group_col}})) +
    scale_color_manual(values = c("red3", "red3"),
                       guide = "none") +
    scale_x_continuous(breaks = seq(1,28,2),
                       expand = expansion(add = 0.25)) +
    scale_y_continuous(labels = ~ifelse(.x <1e3, .x, number(.x / 1e3, 1, suffix = "k")),
                       n.breaks = 10,
                       expand = expansion(mult = c(0, 0.05)),
                       sec.axis = sec_axis(
                         ~rescale(., to = c(0,1)),
                         name = NULL,
                         labels = label_percent()),
                       guide = guide_axis(cap = "none")) +
    coord_cartesian(clip = "off") +
    theme(legend.position.inside = c(0.95,1),
          legend.position = "inside",
          plot.title = element_text(size = 10, hjust = 0),
          axis.title = element_text(size = 9),
          legend.text = element_text(margin = margin(l = 2)),
          legend.justification = c(1, 1),
          legend.key.size = unit(0.75, "lines"),
          legend.key.spacing.y = unit(0.5, "lines"),
          axis.text = element_text(size = 8),
          # plot.margin = margin(t = 2, unit = "lines"),
          axis.line.y.right = element_line(color = "red3"),
          axis.text.y.right = element_text(color = "red3"),
          axis.ticks.y.right = element_line(color = "red3"),
          axis.title.y.right = element_text(color = "red3"),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
    )
}

df_msmr <- df_msmr_tenk10k |> filter(mr == TRUE)

(p_mr_gwas <- mk_plot_tally(magma_gene, p_title = "TenK10K single-cell MR vs. MAGMA"))
(p_mr_eqtlgen <- mk_plot_tally(mr_eqtlgen, 
                              p_title = "TenK10K single-cell MR vs. eQTLgen bulk whole blood MR",
                        label = list(mr_only = "TenK10K only", mr_other = "TenK10K & eQTLgen"))
)

# Example: MAGMA (eczema)
# helper function for plotting

make_plot_locus <- function(
    df_assoc, df_ld, lead_var, lead_var_pos,
    chr, x_ins, start, end,
    col_lead_var = "#6A1B9A",
    col_clump_var = "red3",
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
      # rsq_cat = rev(paletteer_dynamic("cartography::blue.pal", 5) |> as.character()),
      # rsq_cat = rev(c("#E65100", "#F9A825", "#76FF03", "#18FFFF", "#1A237E")),
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
      geom_vline(xintercept = lead_var_pos, color = col_lead_var, linewidth = 0.5,
                 linetype = "longdash") + 
      geom_point(aes(fill = "white"),
                 shape = "circle filled",
                 data = ~.x[variant_id == lead_var], stroke = 1,
                 color = col_lead_var, size = 3.5) +
      geom_point(aes(fill = col_clump_var),
                 shape = "circle filled",
                 data = ~.x[ins == TRUE], alpha = 1,
                 size = 3, stroke = 0) +
    scale_fill_identity(labels = c(label$clump_var, label$lead_var),
                        breaks = c(col_clump_var, "white"),
                         name = NULL) +
    scale_color_viridis_b(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        right = FALSE,
        direction = 1,
        begin = 0.2, end = 1,
        option = "viridis",
        # colours = col_scheme$rsq_cat,
        na.value = col_scheme$rsq_missing,
        name = bquote("LD"~italic(R) ^ 2)
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

trait_name <- "eczema"
cell_type <- "CD14_Mono"
cell_type_print <- "CD14 Mono"
chrNumber <- 12

gene_name_1 <- "ENSG00000139192"
gene_abname_1 <- "TAPBPL"

# to change
gene_name_2 <- "ENSG00000215039"
gene_abname_2 <- "CD27-AS1"
gene_name_3 <- "ENSG00000139190"
gene_abname_3 <- "VAMP1"
gene_name_4 <- "ENSG00000111639"
gene_abname_4 <- "MRPL51"
gene_name_5 <- "ENSG00000010292"
gene_abname_5 <- "NCAPD2"

df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")
# df_msmr_ins <- fread("results/aggregate/tenk10k_phase1.snps")

genes <- c(gene_name_1, gene_name_2, gene_name_3, gene_name_4, gene_name_5)
df_ins <- fread(glue("results/sensitivity/smr/tenk10k_phase1/{cell_type}/{trait_name}/all_chr.snps4msmr.list")) |> 
  filter(gene %in% genes)

df_genes <- df_gene_annot |> 
  filter(ensembl_gene_id %in% genes) |>
  mutate(prefix = glue("results/smr_locus/tenk10k_phase1/{cell_type}/{trait_name}/{ensembl_gene_id}"),
         ld = map(glue("{prefix}.ld"), fread),
         eqtl = map(glue("{prefix}.eqtl.tsv"), fread),
         gwas = map(glue("{prefix}.gwas.tsv"), fread),
         gene_col = paletteer::paletteer_d("ggthemes::Classic_10", length(hgnc_symbol)) |> as.character(),
         ins = map(ensembl_gene_id, ~df_ins[gene == .x]$variant),
         lead_var = map2_chr(eqtl, ins, ~.x[variant_id %in% .y][which.min(p)]$variant_id),
        lead_var_pos = map2_dbl(eqtl, ins, ~.x[variant_id %in% .y][which.min(p)]$pos),
        plot = pmap(list(eqtl, ld, lead_var, lead_var_pos, chr, ins, start-1e5, end+1e5, gene_col, gene_col), make_plot_locus)
        )

# stack the plot and ensure that xlims are shared

xlims <- range(c(df_genes$start - 1e5, df_genes$end + 1e5))

# get gwas data
gwas_prefix <- glue("results/gwas_locus/{trait_name}/tenk10k_phase1/{chrNumber}_{xlims[1]}_{xlims[2]}")
df_gwas <- fread(glue("{gwas_prefix}.gwas.tsv"))
setnames(df_gwas, "SNP", "variant_id")

ld_gwas <- fread(glue("{gwas_prefix}.ld"))
trait_label <- df_trait_map[trait_id == trait_name, label]
p_gwas <- make_plot_locus(
  df_gwas, ld_gwas, 
  lead_var = ld_gwas[1,SNP_A],
  lead_var_pos = ld_gwas[1,BP_A],
  col_lead_var = "black",
  col_clump_var = "black",
  chr = chrNumber,
  x_ins = ld_gwas[1,SNP_A], start = xlims[1], end = xlims[2]
) +
 annotate(
    "label", label = glue("GWAS: {trait_label}"), x = -Inf, y = Inf, vjust = 1, hjust = 0,
    size = 10/.pt, color = "black", fill = "white", linewidth = 0
  )

# get gene tracks using locuszoomr
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
# remap dplyr select & filter
select <- dplyr::select
filter <- dplyr::filter
loc <- locus(xrange = xlims, seqname = chrNumber, ens_db = "EnsDb.Hsapiens.v86")

# add gene colors
gene_colors <- deframe(df_genes |> dplyr::select(hgnc_symbol, gene_col))
loc$TX$gene_col <- gene_colors[loc$TX$gene_name]
loc$TX$exon_col <- gene_colors[loc$TX$gene_name]
loc$TX$exon_border <- gene_colors[loc$TX$gene_name]

p_genetrack <- gg_genetracks(loc, filter_gene_name = df_genes$hgnc_symbol, italics = TRUE)
list_plot <- c(
  list(p_gwas + theme(axis.text.x = element_blank(), legend.position = "none")),
  pmap(list(df_genes$plot, df_genes$hgnc_symbol, df_genes$gene_col, seq_along(df_genes$plot)),
  \(p, gene, col, i) {
    p + ggtext::geom_richtext(
          data = data.frame(x = -Inf, y = Inf),
          aes(x = x, y = y),
          label = glue("{cell_type_print}:<br>*{gene}*"),
          vjust = 1, hjust = 0,
          size = 10 / .pt, color = col,
          fill = "white", label.size = 0
        ) + theme(axis.text.x = element_blank(),
                  legend.position = if (gene == "VAMP1") "inside" else "none")
  }),
  list(p_genetrack + coord_cartesian(clip = "off") + scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
)

(p_locus <- wrap_plots(
  list_plot,
  ncol = 1,
  byrow = TRUE,
  guides = "keep", axis_titles = "collect",
) & 
  theme(panel.grid.major.x = element_line(color = "gray90"),
        panel.ontop = FALSE,
        legend.position.inside = c(1, 0),
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA),
        legend.justification.inside = c(1, 0)) &
  scale_x_continuous(limits = xlims, labels = label_number(0.1, scale = 1e-6, suffix = " Mb"))
)

# P-value MR vs MAGMA vs coloc 
df_plot <- df_msmr_tenk10k |> 
  filter(phenotype == trait_name, biosample == .env$cell_type, probeID %in% genes)
df_plot[df_magma, `:=`(lfdr_magma = i.lfdr, p_magma = i.P),
        on = c("probeID" = "GENE", "phenotype")]
make_plot_comparison <- function(df, y_var, annot_label, threshold = NULL) {
  df <- df |> mutate(.y_var = {{ y_var }},
                     .significant = if (!is.null(threshold)) .y_var >= threshold else TRUE)

  p <- df |>
    ggplot(aes(x = Gene, y = .y_var, fill = Gene, alpha = .significant)) +
    geom_col(show.legend = TRUE) +
    scale_fill_manual(
      values = gene_colors,
      name = NULL,
      guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
    theme_bw() +
    labs(x = NULL, title = annot_label)

  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed", color = "black")
  }
  p
}

p_magma   <- make_plot_comparison(df_plot, -log10(lfdr_magma),       bquote(-log[10] ~ LFDR[MAGMA]), threshold = -log10(0.05))
p_msmr    <- make_plot_comparison(df_plot, -log10(lfdr_msmr_pheno),  bquote(-log[10] ~ LFDR[MR]),    threshold = -log10(0.05))
p_coloc   <- make_plot_comparison(df_plot, coloc_pph4,                bquote(PP[H4] ~ Coloc),         threshold = 0.8)
p_mvcoloc <- make_plot_comparison(df_plot, mvcoloc_pph4,              bquote(PP[H4] ~ "MVColoc"),     threshold = 0.8)

(p_comparison <- wrap_plots(
  list(p_magma, p_msmr, p_coloc, p_mvcoloc),
  nrow = 1,
  byrow = TRUE,
  guides = "collect"
) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 1)),
        legend.text = element_text(size = 9, face = "italic", margin = margin(l = 1)),
        legend.key.spacing.x = unit(0.5, "lines"),
      legend.position = "bottom"))


# combine all plots
plots_example_magma <- wrap_plots(
  list(p_locus, p_comparison),
  ncol = 1,
  byrow = TRUE,
  heights = c(5, 1),
  guides = "keep"
)  +
  plot_annotation(tag_levels = list(c("b", "", "", "", "", "", "", "c"))) &
  theme(plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))


# ggsave("figures/gene_example/magma.png", plots, width = 5, height = 8,
#         device = ragg::agg_png, scaling = 0.8)

# Example: cell-type specificity (BANK1)

trait_name  <- "sle"
trait_label <- "SLE"
chrNumber   <- 4

cell_types <- tribble(
  ~study, ~cell_type,   ~cell_type_print,
  "eqtlgen2020", "bulk_wb", "Whole blood (eQTLgen)",
  "tenk10k_phase1", "B_naive", "B naive",
  "tenk10k_phase1", "B_memory", "B memory",
  "tenk10k_phase1", "CD4_TCM", "CD4 TCM",
  "tenk10k_phase1", "CD8_TEM", "CD8 TEM",
  "tenk10k_phase1", "cDC2", "cDC2"
)

gene_name   <- "ENSG00000153064"
gene_abname <- "BANK1"

df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

# --- load eQTL and LD data per cell type ---

df_ct <- cell_types |>
  mutate(
    prefix  = glue("results/smr_locus/{study}/{cell_type}/{trait_name}/{gene_name}"),
    ld      = map(glue("{prefix}.ld"),       fread),
    eqtl    = map(glue("{prefix}.eqtl.tsv"), fread),
    gene_col = c("gray40", paletteer::paletteer_d("ggthemes::Classic_10", n() -1) |> as.character())
  )

# shared x limits from gene annotation
gene_info <- df_gene_annot |> filter(ensembl_gene_id == gene_name)
xlims <- c(gene_info$start - 1e5, gene_info$end + 1e5)

# instruments (lead eQTL per cell type = top p variant)
df_ct <- df_ct |>
  mutate(
    lead_var     = map_chr(eqtl, ~.x[which.min(p)]$variant_id),
    lead_var_pos = map_dbl(eqtl, ~.x[which.min(p)]$pos),
    plot = pmap(
      list(eqtl, ld, lead_var, lead_var_pos, cell_type_print, gene_col),
      \(eqtl, ld, lead_var, lead_var_pos, ct_print, col) {
        make_plot_locus(
          df_assoc     = eqtl,
          df_ld        = ld,
          lead_var     = lead_var,
          lead_var_pos = lead_var_pos,
          col_lead_var = col,
          col_clump_var = col,
          chr          = chrNumber,
          x_ins        = lead_var,
          start        = xlims[1],
          end          = xlims[2]
        ) +
          ggtext::geom_richtext(
            data = data.frame(x = -Inf, y = Inf),
            aes(x = x, y = y),
            label = glue("{ct_print}:<br>*{gene_abname}*"),
            vjust = 1, hjust = {if (ct_print %in% c("CD4 TCM", "CD8 TEM")) 1 else 0},
            x = {if (ct_print %in% c("CD4 TCM", "CD8 TEM")) Inf else -Inf},
            size = 10 / .pt, color = col,
            fill = "white", label.size = 0
          ) +
          theme(axis.text.x = element_blank(), legend.position = "none") +
          scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
      }
    )
  )

# --- GWAS panel ---

gwas_prefix <- glue("results/gwas_locus/{trait_name}/tenk10k_phase1/{chrNumber}_{xlims[1]}_{xlims[2]}")
df_gwas <- fread(glue("{gwas_prefix}.gwas.tsv"))
setnames(df_gwas, "SNP", "variant_id")
ld_gwas <- fread(glue("{gwas_prefix}.ld"))

p_gwas <- make_plot_locus(
  df_gwas, ld_gwas,
  lead_var      = ld_gwas[1, SNP_A],
  lead_var_pos  = ld_gwas[1, BP_A],
  col_lead_var  = "black",
  col_clump_var = "black",
  chr           = chrNumber,
  x_ins         = ld_gwas[1, SNP_A],
  start         = xlims[1],
  end           = xlims[2]
) +
  annotate(
    "label", label = glue("GWAS: {trait_label}"),
    x = -Inf, y = Inf, vjust = 1, hjust = 0,
    size = 10 / .pt, color = "black", fill = "white", linewidth = 0
  ) +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# --- gene track ---

library(locuszoomr)
library(EnsDb.Hsapiens.v86)

loc <- locus(xrange = xlims, seqname = chrNumber, ens_db = "EnsDb.Hsapiens.v86")

cell_type_colors <- deframe(df_ct |> dplyr::select(cell_type, gene_col))
p_genetrack <- gg_genetracks(loc, filter_gene_name = gene_abname, italics = TRUE)


# comparison of MR effect estimates
df_plot_mr <- df_msmr_tenk10k |> 
  filter(probeID == gene_name, phenotype == trait_name, biosample %in% df_ct$cell_type) |> 
  # add eQTLgen results
  bind_rows(
    df_msmr_eqtlgen |> 
      filter(probeID == gene_name, phenotype == trait_name)
  )

(p_plot_mr <- df_plot_mr |>
  right_join(df_ct |> dplyr::select(cell_type, cell_type_print), by = c("biosample" = "cell_type")) |>
  ggplot(aes(x = factor(cell_type_print, levels = df_ct$cell_type_print), y = -log10(lfdr_msmr_pheno))) +
  geom_col(aes(fill = factor(cell_type_print, levels = df_ct$cell_type_print), alpha = -log10(lfdr_msmr_pheno) >= -log10(0.05)), name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  scale_x_discrete(breaks = df_ct$cell_type_print) +
  scale_fill_manual(values = deframe(df_ct |> select(cell_type_print, gene_col)), name = NULL,
                    breaks = df_ct$cell_type_print,
                    guide = guide_legend(override.aes = list(size = 1))) +
  labs(x = NULL, y = bquote(-log[10]~"LFDR"[MR])) +
  theme(legend.key.size = unit(1, "lines"),
        axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
)

# --- assemble ---
list_plot <- c(
  list(p_gwas),
  df_ct$plot,
  list(p_genetrack + coord_cartesian(clip = "off") +
         scale_y_continuous(expand = expansion(mult = c(0, 0.1))))
)

(p_locus <- wrap_plots(
  list_plot,
  ncol = 1,
  byrow = TRUE,
  guides = "collect",
  axis_titles = "collect"
) &
  theme(panel.grid.major.x = element_line(color = "gray90"),
        panel.ontop = FALSE,
      legend.position = "right") &
  scale_x_continuous(limits = xlims, labels = label_number(0.1, scale = 1e-6, suffix = " Mb"))
)

plots_example_specificity <- (wrap_elements(full = p_locus, clip = FALSE) / p_plot_mr) +
  plot_layout(heights = c(5, 1)) +
  plot_annotation(tag_levels = list(c("e", "f"))) &
  theme(plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))

# Assemble all together
plots_combined <- wrap_plots(
  list(
    wrap_elements(full = p_mr_gwas),
    wrap_elements(full = p_mr_eqtlgen),
    wrap_elements(full = plots_example_magma),
    wrap_elements(full = plots_example_specificity)
  ),
  ncol = 2,
  byrow = TRUE,
  guides = "keep",
  heights = c(1, 5)
) +
  plot_annotation(tag_levels = list(c("a", "d", "", "")))  &
  theme(plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))

ggsave("figures/strict/mr_example/combined.png", plots_combined, width = 8, height = 10,
       device = ragg::agg_png, scaling = 0.8)


# fill in data
# N gene trait not in MAGMA
df_magma[, mr_gene := FALSE]
df_magma[unique(df_msmr_tenk10k[mr == TRUE, .(phenotype, probeID)]),
         mr_gene := TRUE,
         on = c("phenotype", "GENE" = "probeID")]
df_magma[, magma_gene := lfdr < 0.05]

df_msmr_tenk10k[df_magma, magma := i.lfdr < 0.05, on = c("probeID" = "GENE", "phenotype" = "phenotype")]

df_magma_mr <- df_msmr_tenk10k[,
  .(mr = max(mr, na.rm = TRUE),
    magma = max(magma, na.rm = TRUE)),
  by = .(phenotype, probeID)]

tab_magma_mr <- df_magma_mr[,.N, by = .(mr = as.logical(mr), magma = as.logical(magma))]

# number of genes
tab_magma_mr[mr == TRUE & !is.na(magma)] |> 
  mutate(prop = N / sum(N))

# median N genes per trait
tab_magma_mr_by_pheno <- df_magma_mr[mr == 1 & !is.na(magma), .N, by = .(magma, phenotype)] |> 
  mutate(prop = N / sum(N), .by = phenotype)
summary(tab_magma_mr_by_pheno[magma == 0, N])

tab_magma_mr[magma == TRUE & !is.na(mr)] |> 
  mutate(prop = N / sum(N))

# comparison with eQTLgen data
df_msmr_tenk10k[df_msmr_eqtlgen, mr_eqtlgen := i.lfdr_msmr_pheno < 0.05, on = c("probeID", "phenotype")]
df_msmr_tenk10k[is.na(mr_eqtlgen), mr_eqtlgen := FALSE]
df_eqtlgen_tenk10k <- df_msmr_tenk10k[,
  .(tenk10k = max(mr, na.rm = TRUE),
    eqtlgen = max(mr_eqtlgen, na.rm = TRUE)),
  by = .(phenotype, probeID)]

tab_eqtlgen_tenk10k <- df_eqtlgen_tenk10k[,.N, by = .(tenk10k = as.logical(tenk10k), eqtlgen = as.logical(eqtlgen))]

# number of genes
tab_eqtlgen_tenk10k[tenk10k == TRUE & !is.na(eqtlgen)] |> 
  mutate(prop = N / sum(N))

# median N genes per trait
tab_eqtlgen_tenk10k_by_pheno <- df_eqtlgen_tenk10k[tenk10k == 1 & !is.na(eqtlgen), .N, by = .(eqtlgen, phenotype)] |> 
  mutate(prop = N / sum(N), .by = phenotype)
summary(tab_eqtlgen_tenk10k_by_pheno[eqtlgen == 0, N])

tab_eqtlgen_tenk10k[eqtlgen == TRUE & !is.na(tenk10k)] |> 
  mutate(prop = N / sum(N))

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

ggsave("figures/strict/mr_example/unique_celltype_assoc.png", p_celltype_specific, width = 7, height = 4,
       device = ragg::agg_png, scaling = 0.8)
s