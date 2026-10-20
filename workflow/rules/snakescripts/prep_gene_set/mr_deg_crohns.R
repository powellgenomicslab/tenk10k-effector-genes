# preprocessing script with stringent MSMR results as the baseline results

# preprocess main results
library(data.table)
library(tidyverse)
library(arrow)
library(fs)
library(readxl)
library(qvalue)

dir_out <- snakemake@output[[1]]
df_cell_map <- fread("resources/metadata/cell_map.tsv")
df_trait_map_all <- fread("metadata/trait.tsv")
df_trait_map <- filter(df_trait_map_all, include)
df_gene_annot <- fread("resources/misc/gencode.v44.gene_type.tsv")

setDT(df_trait_map)
phenotypes <- df_trait_map$trait_id

cat_order <- fread("metadata/trait_category.tsv") %>%
  pull(cat_order)

# df_msmr_tenk10k <- read_parquet("results/aggregate/tenk10k_phase1.msmr.parquet.gz") %>%
#   filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
#   setDT(key = c("biosample", "phenotype", "probeID"))

df_msmr_tenk10k <- read_parquet("results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz") |>
  filter(!is.na(p_SMR_multi), b_GWAS != 0, b_SMR != 0) |>
  filter(phenotype %in% phenotypes) |>
  mutate(qval_msmr_pheno = qvalue(p_SMR_multi)$qvalues,
         lfdr_msmr_pheno = qvalue(p_SMR_multi)$lfdr,
         pbh_msmr_pheno = p.adjust(p_SMR_multi, "BH"),
         .by = "phenotype") |>
  left_join(df_cell_map %>% select(biosample = wg2_scpred_prediction, cell_type, major_cell_type)) %>%
  inner_join(df_trait_map %>%
    select(
      phenotype = trait_id, pheno_label = label,
      pheno_cat = cat_rev, supercategory
    )) %>%
  group_by(phenotype) %>%
  mutate(
    cell_type = factor(cell_type, df_cell_map$cell_type),
    major_cell_type = factor(major_cell_type, unique(df_cell_map$major_cell_type)),
    pheno_cat = factor(pheno_cat, cat_order)
  ) %>%
  setDT(key = c("biosample", "phenotype", "probeID"))

  gene_universe <- unique(df_msmr_tenk10k$probeID)

  # annotate gene and phenotype
  df_msmr_tenk10k[df_gene_annot, gene_type := i.gene_type, on = c("probeID" = "ensembl_gene_id")]


# Prepare Crohn's DEG for downstream analysis
crohns_dir <- "resources/crohns_case_study"

# manually add the major cell types, then read in the updated file
deg_cell_map <- read.csv(paste0(crohns_dir, "/deg/cell_annotation_or_features/deg_celltype_groups.csv")) %>%
  filter(major_cell_type %in% unique(df_cell_map$major_cell_type))

# Find sig deg
df_deg <- readxl::read_xlsx(paste0(crohns_dir, "/deg/Kongetal2023_supplementary/Kongetal2023_Crohns-sc_DEG.xlsx"), sheet = 1) %>%
  # remove random cols and cont DEG model
  select(1:10, -`Continuous DE coefficients`, -`Continuous DE coefficients p value`, -`Continuous FDR`) %>%
  # get immune cells only
  rename(scRNAseq_cellid = `Cell subset`) %>%
  filter(grepl(pattern = "Immune", scRNAseq_cellid)) %>%
  mutate(scRNAseq_cellid = gsub("Immune.", "", as.character(scRNAseq_cellid)),
         # get the transformed p values
        disc_p_transform = -log10(`Discrete DE coefficients p value`)*sign(`Discrete DE coefficients`),
        # Renaming the contrast and location
        Location = ifelse(Location == "CO", "Colon", "Terminal Ileum") %>% as.factor(),
        Contrast = ifelse(Contrast == "Infl vs. Heal", "Inflamed Tissue vs Healthy", "Non-Inflamed Tissue vs Healthy") %>% as.factor()) %>%
  # annotate with major cell types, removing unmatched entries
  inner_join(deg_cell_map, by = "scRNAseq_cellid") %>%
  # filter to significant only and get rid of cell types not in tenk10k and sign of both discrete and continous are matched
  filter(`Discrete FDR` < 0.05, Location == "Colon")

# Create columns for eviedence criteria and count

evidence_criteria <- list(
  mr = expression(lfdr_msmr_pheno < 0.05)
  # mr_heidi = expression(lfdr_msmr_pheno < 0.05 & p_HEIDI >= 0.05),
  # coloc_single = expression(coloc_pph4 >= 0.8),
  # coloc_multi = expression(mvcoloc_pph4 >= 0.8),
  # mr_plus = expression(lfdr_msmr_pheno < 0.05 & (p_HEIDI >= 0.05 | coloc_single | coloc_multi))
)

for (e in names(evidence_criteria)) {
  df_msmr_tenk10k[, (e) := eval(evidence_criteria[[e]])]
}

df_mr_deg <- df_msmr_tenk10k %>%
  filter(phenotype == "crohns" & mr) %>%
  select(Gene, probeID, biosample, b_SMR, se_SMR, p_SMR, p_SMR_multi, p_HEIDI, major_cell_type) %>%
  inner_join(df_deg, by = c("Gene", "major_cell_type")) %>%
  mutate(concordant_mr_deg = sign(b_SMR) == sign(`Discrete DE coefficients`)) %>%
  mutate(concordance_group = case_when(
            all(concordant_mr_deg) ~ "All concordant",
            all(!concordant_mr_deg) ~ "All discordant",
            any(concordant_mr_deg) & any(!concordant_mr_deg) ~ "Mixed concordance"
          ), .by = c("Gene"))

filters <- list(
  all = expression(TRUE),
  all_no_hla = expression(!grepl("^HLA-", Gene)),
  concordant = expression(concordance_group == "All concordant"),
  concordant_no_hla = expression(!grepl("^HLA-", Gene) & concordance_group == "All concordant"),
  discordant = expression(concordance_group == "All discordant"),
  discordant_no_hla = expression(!grepl("^HLA-", Gene) & concordance_group == "All discordant")
)

for (f in names(filters)) {
  df_out <- df_mr_deg %>% filter(eval(filters[[f]]))
  for (b in c(unique(df_out$biosample), "overall")) {
    if (b == "overall") {
      df <- df_out
    } else {
      df <- df_out[biosample == b]
    }
    if (nrow(df) == 0) {
      next
    }
    gene_set <- df$probeID
    outfile <- file.path(dir_out, b, paste0(f, ".txt"))
    dir_create(dirname(outfile))
    writeLines(gene_set, outfile)
  }
}
