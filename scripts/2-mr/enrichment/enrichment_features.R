# source("scripts/preprocess_strict.R")
library(data.table)
library(tidyverse)
library(ragg)
library(patchwork)
library(arrow)
library(tidygraph)
library(ggraph)
library(paletteer)
library(graphlayouts)
library(tidytext)
library(ggrepel)

df_msmr_tenk10k <- read_parquet("results/preprocessed/tenk10k_phase1.v2.parquet.gz")
df_trait_map <- readxl::read_excel("resources/metadata/trait_metadata_curated.xlsx") |> 
  filter(include) |> 
  setDT()
df_cell_map <- fread("resources/metadata/cell_map.tsv")
cat_order <- readxl::read_xlsx("resources/metadata/trait_metadata_curated.xlsx",
  sheet = "trait_category_order"
) %>%
  pull(cat_order)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")
# count independent phenotypes

INPUT <- list(
  gen_cor = "results/aggregate/tenk10k_phase1.gen_cor.ldak.tsv",
  mr_rrho = "results/rrho/spearman_corr_all_trait_combos_strictmr.tsv"
)
# Create phenotype correlation matrix (same as in calc_indep_clust.R)
# similar analysis using RRHO correlation estimates
df_rrho <- read_tsv(INPUT[["mr_rrho"]]) %>% 
  select(trait1 = p1, trait2 = p2, value = spearman_corr)

# Build square symmetric matrix with 1s on the diagonal
rrho_traits <- union(df_rrho$trait1, df_rrho$trait2)
mat_cor <- bind_rows(
    df_rrho,
    rename(df_rrho, trait1 = trait2, trait2 = trait1),
    tibble(trait1 = rrho_traits, trait2 = rrho_traits, value = 1)
  ) %>%
  pivot_wider(names_from = trait2, values_from = value) %>%
  column_to_rownames("trait1") %>%
  data.matrix() %>%
  .[rrho_traits, rrho_traits]

# Apply to your data (assuming df_msmr_tenk10k is loaded)
# Choose your evidence criteria (e.g., "mr_plus")
gene_summary <- df_msmr_tenk10k[mr == TRUE, .(
  n_phenotypes = length(unique(pheno_label)),
  n_pheno_cat = length(unique(pheno_cat)),
  n_celltypes = length(unique(biosample)),
  phenotypes = list(unique(pheno_label)),
  celltypes = list(unique(biosample))
), by = probeID]

# Function using eigenvalue decomposition
count_effective_tests <- function(pheno_vec, cor_matrix) {
  if (length(pheno_vec) <= 1) return(length(pheno_vec))
  
  pheno_in_matrix <- intersect(pheno_vec, rownames(cor_matrix))
  if (length(pheno_in_matrix) <= 1) return(length(pheno_in_matrix))
  
  sub_cor <- cor_matrix[pheno_in_matrix, pheno_in_matrix]
  
  eigenvalues <- eigen(sub_cor, symmetric = TRUE, only.values = TRUE)$values
  
  # Li and Ji (2005) method
  f <- function(x) {
    if (x <= 0) return(0)
    indicator <- as.numeric(x >= 1)
    fractional <- x - floor(x)
    return(indicator + fractional)
  }
  
  meff <- sum(sapply(eigenvalues, f))
  return(meff)
}

# Apply to gene summary
gene_summary[, n_effective_phenotypes := sapply(phenotypes, function(p) {
  count_effective_tests(p, mat_cor)
})]


enrich <- "string_rank_api"
biosample <- "CD4_CTL"
# enrichment_file <- glue::glue("results/aggregate/enrichment/tenk10k_phase1_mr_strict.mr_strict.{enrich}.tsv.gz")
enrichment_file <- glue::glue("results/enrichment/tenk10k_phase1_mr_strict/mr_strict/{biosample}.{enrich}.tsv")
df_enrich <- fread(enrichment_file)

df_enrich[df_trait_map, `:=`(
  pheno_label = i.label, pheno_cat = i.cat_rev, supercategory = i.supercategory), on = c("phenotype" = "trait_id")]

# annotate terms that contains HLA genes
df_enrich[, l_intersection_genes := str_split(intersection_genes, ",")]
df_enrich[,
  prop_hla := sapply(l_intersection_genes, function(genes) {
                    mean(str_detect(genes, "HLA"))
})]
df_enrich[, term_name := str_replace_all(term_name, "Mixed, incl\\. ", "Mixed: ")]
df_enrich[, maj_hla := prop_hla > 0.5]

sources <- c("GO Process", "KEGG", "STRING clusters")

sources <- c("GO Process", "KEGG", "STRING clusters", "Reactome")
df_enrich |>
  filter(source != "Publications", supercategory == "disease", source %in% sources) |>
  filter(!maj_hla) |>
  # filter(source %in% sources) |>
  ggplot(aes(x = -log10(fdr), color = pheno_cat, y = reorder_within(term_name |> str_trunc(40), -log10(fdr), source))) +
  geom_point() +
  scale_y_reordered() +
  theme_bw() +
  facet_wrap(vars(source), scales = "free_y", space ="free_y", strip.position = "left") +
  theme(legend.position = "bottom", strip.clip = "off")

# get numbers
# N tested (disease only)
df_msmr_tenk10k[supercategory == "disease" & mr == TRUE & biosample == get("biosample", envir = parent.env(environment())), uniqueN(phenotype)]

# N signif
df_enrich[supercategory == "disease" & source %in% sources, .N]

# N terms excluding hla
df_enrich[supercategory == "disease" & source %in% sources & maj_hla == FALSE, .(n_phenotypes = uniqueN(phenotype)), by = .(term_name, term_id)][order(-n_phenotypes)]

# terms with 1 phenotype only
df_enrich[supercategory == "disease" & source %in% sources & maj_hla == FALSE] |> 
  group_by(term_name, term_id) |>
  filter(n() == 1) |> 
  View()


# --- Bipartite network: terms (left) <-> phenotypes (right) ---
df_bip_edges <- df_enrich |>
  as_tibble() |>
  filter(source %in% sources, supercategory == "disease", !maj_hla) |>
  mutate(
    term_label = str_trunc(term_name, 40),
    log10_fdr  = -log10(fdr)
  ) |>
  select(term_label, pheno_label, pheno_cat, source, log10_fdr)

# Flag phenotype-specific terms (connected to exactly one phenotype)
pheno_specific_terms <- df_bip_edges |>
  distinct(term_label, pheno_label) |>
  count(term_label) |>
  filter(n == 1) |>
  pull(term_label)

df_bip_edges <- df_bip_edges |>
  mutate(is_specific = term_label %in% pheno_specific_terms)

# Shared pheno_cat palette; non-specific gets "grey70"
pheno_cats  <- sort(unique(df_bip_edges$pheno_cat))
# cat_palette <- setNames(
#   paletteer_d("ggthemes::Tableau_10", n = length(pheno_cats)),
#   pheno_cats
# )

# Sort terms: group by source, specific before non-specific within each source
# One row per term name — for terms in multiple pheno_cats pick the specific one
term_nodes <- df_bip_edges |>
  arrange(source, desc(is_specific), pheno_cat, term_label) |>
  distinct(name = term_label, .keep_all = TRUE) |>
  select(name, source, pheno_cat, is_specific) |>
  mutate(node_type = "term")

pheno_nodes <- df_bip_edges |>
  distinct(name = pheno_label, pheno_cat) |>
  arrange(pheno_cat, name) |>
  mutate(node_type = "phenotype", source = NA_character_, is_specific = TRUE)

bip_nodes <- bind_rows(term_nodes, pheno_nodes)

# Store intended term order before tbl_graph reorders nodes
term_order <- setNames(seq_along(term_nodes$name), term_nodes$name)

bip_edges <- df_bip_edges |>
  select(from = term_label, to = pheno_label, pheno_cat, source, log10_fdr, is_specific)

g_bip <- tbl_graph(nodes = bip_nodes, edges = bip_edges, directed = FALSE) |>
  mutate(type = node_type == "phenotype")

# Manual layout
layout_bip <- g_bip |> create_layout(layout = "bipartite")
layout_bip$x <- ifelse(layout_bip$node_type == "term", 0, 1)

term_idx  <- which(layout_bip$node_type == "term")
pheno_idx <- which(layout_bip$node_type == "phenotype")

# Assign y: uniform inter-node spacing across all terms, with a fixed gap
# inserted at each source boundary so brackets never overlap
n_terms     <- nrow(term_nodes)
n_sources   <- length(unique(term_nodes$source))
gap_steps   <- 2L  # gap between source groups = equivalent of this many node spacings
total_steps <- (n_terms - 1) + gap_steps * (n_sources - 1)
step        <- 1 / total_steps

term_layout_info <- term_nodes |>
  mutate(rank_within = row_number()) |>
  group_by(source) |>
  mutate(source_idx = cur_group_id(), pos_within = row_number() - 1L) |>
  ungroup() |>
  mutate(
    # cumulative position: node rank + gaps accumulated before this source group
    cum_pos = (rank_within - 1L) + gap_steps * (source_idx - 1L),
    y_val   = 1 - cum_pos * step
  )

term_name_to_y <- setNames(term_layout_info$y_val, term_layout_info$name)
layout_bip$y[term_idx]  <- term_name_to_y[layout_bip$name[term_idx]]
layout_bip$y[pheno_idx] <- seq(0, 1, length.out = length(pheno_idx))

# Source bracket annotations: span ALL terms (specific + non-specific) per source
source_labels <- layout_bip[term_idx, ] |>
  group_by(source) |>
  summarise(y_mid = mean(y), y_min = min(y), y_max = max(y), .groups = "drop")

layout_bip_edges <- attr(layout_bip, "graph") |>
  activate(edges) |>
  as_tibble()

edges_nonspecific <- layout_bip_edges |> filter(!is_specific)
edges_specific    <- layout_bip_edges |> filter(is_specific)

cat_palette <- paletteer::paletteer_d("ggthemes::Tableau_10", n = length(cat_order)) |> 
  as.character() |> set_names(cat_order)
  
p_bip <- ggraph(layout_bip) +
  geom_edge_bend(
    aes(edge_width = log10_fdr, filter = !is_specific),
    colour = "grey90", alpha = 0.4, strength = 0.2, show.legend = FALSE
  ) +
  geom_edge_bend(
    aes(edge_width = log10_fdr, colour = pheno_cat, filter = is_specific),
    alpha = 0.7, strength = 0.2, show.legend = FALSE
  ) +
  # Non-specific term nodes (background)
  geom_node_point(
    data = \(d) filter(d, node_type == "term", !is_specific),
    shape = 22, size = 2.5, fill = "grey80", color = "grey60", stroke = 0.3
  ) +
  # Specific term nodes (coloured squares)
  geom_node_point(
    data = \(d) filter(d, node_type == "term", is_specific),
    aes(fill = pheno_cat),
    shape = 22, size = 3, color = "grey20", stroke = 0.3
  ) +
  # Phenotype nodes (coloured circles)
  geom_node_point(
    data = \(d) filter(d, node_type == "phenotype"),
    aes(fill = pheno_cat),
    shape = 21, size = 3, color = "grey20", stroke = 0.3
  ) +
  # Add an invisible point layer just for the fill legend
  geom_point(
    data = \(d) filter(d, node_type == "phenotype"),
    aes(x = x, y = y, fill = pheno_cat),
    shape = 22, size = 0, colour = NA
  ) +
  scale_fill_manual(
    name = NULL,
    values = cat_palette,
    guide = guide_legend(override.aes = list(size = 4, shape = 22, colour = "grey20"))
  ) +
  scale_edge_colour_manual(values = cat_palette, guide = "none") +
  # Non-specific term labels (grey)
  geom_node_text(
    data = \(d) filter(d, node_type == "term", !is_specific),
    aes(label = name),
    hjust = 1, nudge_x = -0.02, size = 2.0, colour = "grey70"
  ) +
  # Specific term labels (black)
  geom_node_text(
    data = \(d) filter(d, node_type == "term", is_specific),
    aes(label = name),
    hjust = 1, nudge_x = -0.02, size = 2.3
  ) +
  geom_node_text(
    data = \(d) filter(d, node_type == "phenotype"),
    aes(label = name),
    hjust = 0, nudge_x = 0.02, size = 2.3
  ) +
  # Source brackets (left margin, specific terms only)
  geom_segment(
    data = source_labels,
    aes(x = -0.5, xend = -0.5, y = y_min, yend = y_max),
    linewidth = 0.5, colour = "grey50", inherit.aes = FALSE
  ) +
  geom_text(
    data = source_labels,
    aes(x = -0.55, y = y_mid, label = source),
    angle = 90, hjust = 0.5, size = 2.8, colour = "grey30", inherit.aes = FALSE
  ) +
  scale_edge_width_continuous(range = c(0.5, 1.5), guide = "none") +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right",
        legend.title.position = "top")

ggsave(
  "figures/enrichment/bipartite_string_enrichment_CD4_CTL.png",
  p_bip, width = 8, height = 8, device = ragg::agg_png, scaling = 1.1
)


# --- PPI network visualisation ---

df_string <- fread("resources/enrichment/string_id_map.tsv")
df_network <- fread(glue::glue("results/enrichment/tenk10k_phase1_mr_strict/mr_strict/{biosample}.{enrich}.ppi_network.tsv"))
df_bipartite <- fread(glue::glue("results/enrichment/tenk10k_phase1_mr_strict/mr_strict/{biosample}.{enrich}.bipartite.tsv"))

# remove housekeeping genes
housekeeping_genes <- fread("resources/misc/HOUNKPE_HOUSEKEEPING_GENES.v2026.1.Hs.tsv") |> 
  filter(STANDARD_NAME == "GENE_SYMBOLS") |> 
  pull(2) |> 
  str_split(",") |> 
  _[[1]]

# Node annotations: number of distinct phenotypes each gene is associated with
gene_summary[df_gene_annot, gene := i.hgnc_symbol, on = c(probeID = "ensembl_gene_id")]
node_annot <- df_bipartite[,
  .(n_phenotypes = uniqueN(phenotype), min_fdr = min(fdr)),
  by = gene
] |> 
  left_join(gene_summary[, .(n_effective_phenotypes, gene)])

# Build node table from all genes appearing in the network
all_genes <- union(df_network$from, df_network$to) |> 
  # exclude housekeeping_genes
  setdiff(housekeeping_genes)
nodes <- tibble(name = all_genes) |>
  left_join(node_annot, by = c("name" = "gene")) |>
  mutate(
    n_phenotypes = replace_na(n_effective_phenotypes, 0L),
    in_bipartite = n_phenotypes > 0
  )

g_ppi <- tbl_graph(
  nodes = nodes,
  edges = df_network |> select(from, to, score),
  directed = FALSE
) |>
  mutate(
    degree         = centrality_degree(weights = score),
    betweenness    = centrality_betweenness(weights = score),
    eigenvector    = centrality_eigen(weights = score),
    pagerank       = centrality_pagerank(weights = score),
    closeness      = centrality_closeness(weights = score)
  )

# Compute layout on the full graph
set.seed(8)

# g_ppi_main <- g_ppi
g_ppi_main <- g_ppi |> filter(degree > 10)
# layout_full <- create_layout(g_ppi_main, layout = "fr", weights = score)
layout_full <- create_layout(g_ppi_main, layout = "backbone", keep = 0.8)

# Extract high-pleiotropy nodes for overlay
pleio_threshold <- quantile(layout_full$n_phenotypes, 0.9)
pleio_nodes <- layout_full |>
  filter(n_phenotypes >= pleio_threshold)

betweenness_threshold <- quantile(layout_full$betweenness, 0.9)
pleio_nodes <- layout_full |>
  filter(n_phenotypes >= 1, betweenness >= betweenness_threshold)

# Edges among high-pleiotropy nodes
edges_full <- g_ppi_main |> activate(edges) |> as_tibble()
pleio_edges <- edges_full |>
  filter(
    layout_full$name[from] %in% pleio_nodes$name &
    layout_full$name[to]   %in% pleio_nodes$name
  ) |>
  mutate(
    x    = layout_full$x[from], y    = layout_full$y[from],
    xend = layout_full$x[to],   yend = layout_full$y[to]
  )

p_graph <- ggraph(layout_full) +
  # Background: full network
  geom_edge_link(alpha = 0.03, colour = "grey80") +
  geom_node_point(size = 0.5, colour = "grey85") +
  # Foreground: edges among high-pleiotropy nodes
  geom_segment(
    data = pleio_edges,
    aes(x = x, y = y, xend = xend, yend = yend, alpha = score),
    colour = "blue3", show.legend = FALSE
  ) +
  # Foreground: high-pleiotropy nodes
  geom_point(
    data = pleio_nodes,
    aes(x = x, y = y, fill = n_phenotypes, size = degree),
    shape = 21, stroke = 0.4, colour = "black"
  ) +
  geom_label_repel(
    data = pleio_nodes,
    aes(x = x, y = y, label = name),
    size = 2.5, max.overlaps = 40, fontface = "italic",
    fill = alpha("white", 0.5), label.size = 0, label.padding = 0.1,
    segment.colour = "grey60", segment.size = 0.2,
    show.legend = FALSE
  ) +
  scale_fill_viridis_c(name = "N independent\nphenotypes", option = "plasma", direction = -1) +
  scale_size_continuous(name = "Weighted degree", range = c(1, 5), guide = "none") +
  scale_alpha_continuous(range = c(0.1, 0.6), guide = "none") +
  theme_graph(base_family = "sans") +
  labs(title = "STRING PPI network (CD4 CTL genes) — high-pleiotropy genes")

ggsave(
  "figures/enrichment/string_high_pleiotropy_network_v2.png",
  p_graph, width = 8, height = 6, device = ragg::agg_png
)

degree_df <- layout_full |>
  as_tibble() |>
  transmute(
    name,
    betweenness,
    n_phenotypes,
    is_mr = n_phenotypes > 0
  )

t_res <- t.test(betweenness ~ is_mr, data = degree_df)

top_genes <- degree_df |>
  filter(is_mr) |>
  slice_max(betweenness, n = 10)

p_violin <- ggplot(degree_df, aes(x = is_mr, y = betweenness)) +
  geom_violin(alpha = 0.5, colour = NA, fill = "grey70",
              data = ~filter(.x, !is_mr)) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, data = ~filter(.x, !is_mr),
               colour = "black", fill = "grey50") +
  geom_violin(alpha = 0.5, colour = NA, fill = "steelblue3",
              data = ~filter(.x, is_mr)) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, data = ~filter(.x, is_mr),
               colour = "black", fill = "steelblue4") +
  geom_point(
    data = top_genes,
    aes(x = is_mr, y = betweenness, fill = n_phenotypes),
    shape = 21, colour = "grey20", stroke = 0.4, show.legend = TRUE
  ) +
  geom_label_repel(
    data = top_genes,
    aes(x = is_mr, y = betweenness, label = name),
    size = 2.5, max.overlaps = 20, fontface = "italic",
    fill = alpha("white", 0.6), label.size = 0,
    segment.colour = "grey50", segment.size = 0.1,
    inherit.aes = FALSE,
    direction = "y", hjust = 1, xlim = c(-Inf, 1.8)
  ) +
  scale_x_discrete(labels = c("FALSE" = "Non-MR genes", "TRUE" = "MR genes")) +
  scale_fill_viridis_c(name = "N independent\nphenotypes", option = "plasma", direction = -1) +
  guides(
    fill = guide_colorbar(title.position = "top")
  ) +
  labs(
    x = NULL,
    y = "Betweenness centrality",
    subtitle = glue::glue(
      "t = {round(t_res$statistic, 2)}, p = {signif(t_res$p.value, 3)} (Welch two-sample t-test)"
    )
  ) +
  theme_classic() +
  theme(legend.position = "right",
        plot.subtitle = element_text(hjust = 0.5))

ggsave(
  "figures/enrichment/string_betweenness_by_mr.png",
  p_violin, width = 6, height = 5, device = ragg::agg_png
)

# --- STRING cluster-focused network ---
# highlight MS-related genes
df_ms <- df_enrich |> 
  filter(source %in% sources, phenotype == "ms")

# highlight_terms <- c("PMID:22190364", "CL:16682", "CL:16585", "hsa04064")
highlight_terms <- unique(df_ms$term_id)

# Compute layout — stress preserves distances well for small highlighted subgraphs
# set.seed(123)
# layout_full <- create_layout(g_ppi, layout = "backbone", keep = 0.1)

# MS MR-significant genes in CD4_CTL
ms_mr_genes <- df_msmr_tenk10k[
  phenotype == "ms" & mr == TRUE & biosample == "CD4_CTL",
  unique(Gene)
]
ms_mr_gene_set <- ms_mr_genes[ms_mr_genes %in% layout_full$name]

# Which of those are in the highlight terms?
ms_term_genes <- df_enrich[
  phenotype == "ms" & term_id %in% highlight_terms,
  .(gene = unlist(str_split(intersection_genes, ","))),
  by = .(term_id, term_name)
] |> unique()

ms_term_genes[, term_label := str_wrap(term_name, 100)]

# Assign highlight term (smallest term for genes in multiple)
ms_term_genes[, term_size := .N, by = term_id]
ms_gene_assign <- ms_term_genes[, .SD[which.min(term_size)], by = gene]

# Build layout data: MR genes + non-MR term genes
# All genes from the enrichment terms (includes non-MR genes)
all_term_genes <- unique(ms_gene_assign$gene)
# Combined gene set: MR genes + term genes
ms_combined_set <- union(ms_mr_gene_set, all_term_genes[all_term_genes %in% layout_full$name])

layout_ms <- layout_full[layout_full$name %in% ms_combined_set, ]
layout_ms$term <- ms_gene_assign$term_label[match(layout_ms$name, ms_gene_assign$gene)]
layout_ms$is_mr <- layout_ms$name %in% ms_mr_gene_set
layout_ms$in_term <- !is.na(layout_ms$term)

# Edges among combined gene set
ms_edges <- edges_full |>
  filter(
    layout_full$name[from] %in% ms_combined_set &
      layout_full$name[to] %in% ms_combined_set
  ) |>
  mutate(x = layout_full$x[from], y = layout_full$y[from],
         xend = layout_full$x[to], yend = layout_full$y[to])

# Palette: enrichment terms + MR-only + term-only categories
ms_term_levels <- sort(unique(na.omit(layout_ms$term)))
n_ms_terms <- length(ms_term_levels)
ms_pal <- c(
  paletteer_d("ggthemes::Tableau_10")[seq_len(n_ms_terms)] |> set_names(ms_term_levels),
  c("Other MR genes" = "grey50")
)

layout_ms$term_fill <- factor(
  case_when(
    !is.na(layout_ms$term) ~ layout_ms$term,
    layout_ms$is_mr ~ "Other MR genes"
  ),
  levels = c(ms_term_levels, "Other MR genes")
)

# Shape: filled circle for MR genes, open circle for non-MR term genes
layout_ms$node_shape <- ifelse(layout_ms$is_mr, "circle filled", "circle open")

p_ms <- ggraph(layout_full) +
  geom_edge_link(alpha = 0.03, colour = "grey80") +
  geom_node_point(size = 0.5, colour = "grey85") +
  geom_segment(
    data = ms_edges,
    aes(x = x, y = y, xend = xend, yend = yend, alpha = score),
    colour = "black", show.legend = FALSE
  ) +
  # Non-MR term genes (open circles, smaller)
  geom_point(
    data = layout_ms[!layout_ms$is_mr, ],
    aes(x = x, y = y, fill = term_fill),
    shape = 21, size = 2.5, alpha = 1, stroke = 0.5
  ) +
  # MR genes (diamonds, larger)
  geom_point(
    data = layout_ms[layout_ms$is_mr, ],
    aes(x = x, y = y, fill = term_fill, shape = "MR genes"),
    size = 4, alpha = 1, stroke = 0.5
  ) +
  geom_label_repel(
    data = layout_ms,
    aes(x = x, y = y, label = name, color = term_fill),
    size = 2.5, max.overlaps = 40, fontface = "plain",
    fill = alpha("white", 0.5),
    show.legend = FALSE, segment.colour = "grey60",
    segment.size = 0.2, label.size = 0
  ) +
  scale_fill_manual(
    name = "Enrichment term", values = ms_pal,
    # labels = ~str_wrap(., 50),
    aesthetics = c("color", "fill"),
  guide = guide_legend(
      ncol = 1,
      theme = theme(legend.title.position = "top"),
      override.aes = list(size = 4, shape = 22)
    )) +
  scale_shape_manual(name = NULL, values = c("MR genes" = "diamond filled"), guide = guide_legend(override.aes = list(fill = "white", stroke = 1.5))) +
  scale_alpha_continuous(range = c(0.1, 0.6), guide = "none") +
  theme_graph(base_family = "sans") +
  # labs(
  #   title = "STRING PPI network (CD4 CTL genes) — MS enriched terms"
  # ) +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.key.spacing.y = unit(0.1, "lines"))

ggsave(
  "figures/enrichment/string_ms_highlight_network.png",
  p_ms, width = 10, height = 8, device = ragg::agg_png,
)

# plot assembly
plots <- wrap_plots(
  list(
    p_bip + expand_limits(x = c(1.5, 1)) +
      theme(legend.box.margin = margin(l = -7, unit = "lines")),
    p_ms),
  ncol = 1,
  heights = c(1.6, 1),
  guide = "keep"
) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"),
        plot.margin = margin(),
      plot.tag.position = c(0.02, 0.97))
  
ggsave("figures/enrichment/string_enrichment_combined.png",
  plots, width = 9, height = 12, device = ragg::agg_png, scaling = 1
)

# get numbers
df_enrich



# Extract non-MHC STRING cluster gene memberships
cluster_genes <- df_enrich[
  source == "STRING clusters" & !str_detect(term_name, "MHC"),
  .(gene = unlist(str_split(intersection_genes, ","))),
  by = .(term_id, term_name)
] |> unique()

# For genes in multiple clusters, pick the most specific (smallest) cluster
cluster_genes[, cluster_size := .N, by = term_id]
cluster_assign <- cluster_genes[, .SD[which.min(cluster_size)], by = gene]

# Shorten cluster labels for plotting
cluster_assign[, cluster_label := str_trunc(
  str_remove(term_name, "^Mixed, incl\\.\\s*"), 50
)]

# Subgraph: only genes that belong to a non-MHC STRING cluster
cluster_gene_set <- unique(cluster_assign$gene)

g_ppi_clusters <- g_ppi |>
  filter(name %in% cluster_gene_set) |>
  mutate(cluster = cluster_assign$cluster_label[match(name, cluster_assign$gene)]) |>
  filter(group_components() == 1)

set.seed(42)

# Compute layout on the full main-component graph so cluster nodes sit in context
layout_full <- create_layout(g_ppi_main, layout = "backbone", keep = 0.2)

# Extract cluster node/edge data for overlay
cluster_node_idx <- which(layout_full$name %in% cluster_gene_set)
layout_cluster <- layout_full[cluster_node_idx, ]
layout_cluster$cluster <- cluster_assign$cluster_label[match(layout_cluster$name, cluster_assign$gene)]

# Edges among cluster genes
edges_full <- g_ppi_main |> activate(edges) |> as_tibble()
cluster_edges <- edges_full |>
  filter(
    layout_full$name[from] %in% cluster_gene_set &
    layout_full$name[to] %in% cluster_gene_set
  ) |>
  mutate(x = layout_full$x[from], y = layout_full$y[from],
         xend = layout_full$x[to], yend = layout_full$y[to])

# Number of unique clusters for palette
n_clusters <- length(unique(layout_cluster$cluster))
cluster_pal <- paletteer_d("ggthemes::Tableau_20")[seq_len(n_clusters)] |>
  set_names(sort(unique(layout_cluster$cluster)))

p_graph <- ggraph(layout_full) +
  # Background: full network edges + nodes in grey
  geom_edge_link(alpha = 0.03, colour = "grey80") +
  geom_node_point(size = 0.5, colour = "grey85") +
  # Foreground: cluster edges
  geom_segment(
    data = cluster_edges,
    aes(x = x, y = y, xend = xend, yend = yend, alpha = score),
    colour = "black", show.legend = FALSE
  ) +
  # Foreground: cluster nodes
  geom_point(
    data = layout_cluster,
    aes(x = x, y = y, size = n_phenotypes, fill = cluster),
    alpha = 1, shape = "circle filled"
  ) +
  geom_label_repel(
    data = layout_cluster,
    aes(x = x, y = y, label = name), colour = "black",
    size = 2.5, max.overlaps = 40, fontface = "italic", fill = alpha("white", 0.5),
    show.legend = FALSE, segment.colour = "grey60", segment.size = 0.2, label.size = 0
  ) +
  scale_fill_manual(name = "STRING cluster", values = cluster_pal) +
  scale_size_continuous(name = "N Phenotypes", range = c(2, 8)) +
  scale_alpha_continuous(range = c(0.1, 0.6), guide = "none") +
  theme_graph(base_family = "sans") +
  guides(fill = guide_legend(ncol = 3,
    theme = theme(legend.title.position = "top"),
    override.aes = list(size = 4, shape = "square filled"))) +
  labs(
    title = "STRING PPI network — non-MHC enriched clusters (CD4 CTL)",
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    # legend.justification = c(1, 0),
    # legend.background = element_blank(),
    # panel.grid = element_blank()
  )

ggsave(
  "figures/enrichment/string_cluster_network.png",
  p_graph, width = 10, height = 10, device = ragg::agg_png
)
