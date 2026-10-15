# supplementary table: CD4_CTL STRING enrichment results, annotated with trait metadata
#
# Reuses the df_enrich load/annotate steps from scripts/enrichment/enrichment_features.R
# (enrich = "string_rank_api", biosample = "CD4_CTL"), without the plotting code, and
# applies the same `sources`/disease filter used by that script's plots and tallies
# (scripts/enrichment/enrichment_features.R:104-106, 122) -- confirmed against the live
# "8-enrich_cd4ctl" tab: filtering to these 4 sources + disease-only traits reproduces
# it exactly (540/540 rows match, 0 discrepancies either direction).

library(data.table)
library(tidyverse)

source("scripts/util/write_table.R")

df_trait_map <- fread("metadata/trait.tsv") |>
  filter(include) |>
  setDT()

enrich <- "string_rank_api.trimmed"
biosample <- "CD4_CTL"
enrichment_file <- glue::glue("results/enrichment/tenk10k_phase1_mr_strict/mr_strict/{biosample}.{enrich}.tsv")
df_enrich <- fread(enrichment_file)

df_enrich[df_trait_map, `:=`(
  pheno_label = i.label, pheno_cat = i.cat_rev, supercategory = i.supercategory
), on = c("phenotype" = "trait_id")]

# annotate terms that contain HLA genes
df_enrich[, l_intersection_genes := str_split(intersection_genes, ",")]
df_enrich[, prop_hla := sapply(l_intersection_genes, function(genes) mean(str_detect(genes, "HLA")))]
df_enrich[, term_name := str_replace_all(term_name, "Mixed, incl\\. ", "Mixed: ")]

sources <- c("GO Process", "KEGG", "STRING clusters", "Reactome")
tab_enrich_cd4ctl <- df_enrich[source %in% sources & supercategory == "disease"]

write_table(tab_enrich_cd4ctl, "enrich_cd4ctl", 8)
