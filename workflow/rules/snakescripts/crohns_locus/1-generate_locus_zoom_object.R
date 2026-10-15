# Purpose: Prepare Locus zoom of eQTL for TenK10K
# Extended Data Figure 9 (GPX1): rule crohns_gpx1_locus_data in workflow/rules/misc.smk
# Requires individual-level genotypes (PLINK) to compute LD; plotted by scripts/8-crohns/plot_figures/2-plot_locus_zoom_fig.R.

# Load renv with locuszoomr and tidyverse
options(renv.config.synchronized.check = FALSE)
# renv::load(project = <renv project with locuszoomr and tidyverse>)

args <- commandArgs(trailingOnly = TRUE)
# Print all arguments to check
print(args)
target_gene <- args[1]
chr <- as.numeric(args[2]) # Convert to number if needed
gene_ens <- args[3]

# target_gene <- "RUNX3"
# chr <- 1 # Convert to number if needed
# gene_ens <- "ENSG00000020633"

cat("Run locus_plot() for :", target_gene, "\n")
cat("on chromosome :", chr, "\n")
####################### Load packages and connect db for locuszoomr functions #######################
library(locuszoomr)
library(tidyverse)
library(data.table)
library(ragg)
#library(cowplot)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
#library(LDlinkR)

ah <- AnnotationHub(localHub = TRUE)
query(ah, c("EnsDb", "Homo sapiens"))
ensDb_v106 <- ah[["AH100643"]]
cat("Load libraries")
####################### Paths #######################
print("Get eQTL for: ")
target_gene
gene_ens
#cellType <- c("BIN", "BMem", "CD4NC", "CD4ET","CD4SOX4", "CD8ET","CD8NC", "CD8S100B","DC","MonoNC", "MonoC", "NK","NKR","Plasma")
cellType <- readLines("resources/crohns_case_study/figures/tenk10k_celltypes.txt")
#cellType <- c("CD4_Naive", "CD4_TCM", "cDC2", "MAIT")
#cellType <- c("CD4_TCM", "cDC2")

gwas_name <- "crohns"
gwas_fp <- paste0("resources/pipeline_ma/", gwas_name, ".ma")
eqtl_folders <- "resources/matrix_eqtl/tenk10k_phase1"
genotype_files="resources/genotypes/tenk10k_phase1"
plink_out = paste0("resources/crohns_case_study/locus_zoom_data/", target_gene)
print("Set paths")

if (!dir.exists(plink_out)) {
  dir.create(plink_out, recursive = TRUE)
}

eqtl_fp <- paste0(eqtl_folders, "/", cellType, "/", "chr", chr, "_meqtl.tsv")

# Only keep cell types that have an eQTL file
eqtl_exists <- file.exists(eqtl_fp)
cat("Skipping cell types with no eQTL file: ", cellType[!eqtl_exists], "\n")
cellType <- cellType[eqtl_exists]
eqtl_fp <- eqtl_fp[eqtl_exists]

print("Set paths")
cat("File paths: ", gwas_fp, "\n", eqtl_folders, "\n", eqtl_fp, "\n", genotype_files, "\n", plink_out)

####################### Read in GWAS data for locuszoomr #######################
# Note you need CHR, BP columns 
gwas <- fread(gwas_fp) 

# need to add the chr, bp columns back for locus_plot(). fyi this is suuuuper slow so maybe use a datatable version?
gwas <- gwas %>%
  tidyr::separate(SNP, into = c("chrom", "pos"), sep = ":", extra = "drop", remove = FALSE) %>%
  dplyr::rename(MarkerID = SNP) # So it matches the eQTL results 

print("Read in GWAS file & rename SNP")
####################### Set the gene symbol and ensembl ID, chromsome and (later) cell type of interest #######################
# Read in SAIGE or Matrix eQTL files for each cell type, only for the pre-set chromosome of interest

locus_input <- lapply(eqtl_fp, fread)
# Create list of the eQTL files across all cell types, and ensure to filter for your gene of interest

# for SAIGE eQTL
# locus_input <- locus_input %>%
#   map(~ .x %>% dplyr::rename(p = "p.value") %>%
#   dplyr::filter(gene == gene_ens))
print("Read in eQTL files")
# for Matrix eQTL: 
locus_input_filtered <- locus_input %>%
  map(~ .x %>% tidyr::separate(SNP, into = c("chrom", "pos"), sep = ":", extra = "drop", remove = FALSE) %>%
  dplyr::rename(p = "p-value", MarkerID = SNP) %>%
  dplyr::filter(gene == gene_ens))


print("Process eQTL files for locus input")
# name the list of eQTL files 
names(locus_input_filtered) <- cellType

print("Set names")
print("Test")

head(locus_input_filtered[[1]])

# add the GWAS to the eQTL file list 
locus_input_filtered[["GWAS"]] <- gwas
print("Add GWAS locus to list")

# Run locus() function on everything in locus_input (NAMED LIST of eQTLs and GWASes, where eQTL SNPs filtered to the signals for a given eGene)
loclist <- lapply(locus_input_filtered, locus, gene = target_gene, ens_db = "EnsDb.Hsapiens.v86", flank = 5e5)
print("Run locus() function on eQTL and GWAS. Names of the list: ")
names(loclist)

# save(gwas_loc, loclist, gene, file = "resources/crohns_case_study/figures/locus_zoom_plot_objects.RData")

print("Run PLINK")

############### GET PLINK DATA 
#plink_chr15_renamed_201010.bim
# Loop through every item in your locus plot data list
for (cell in names(loclist)) {
    # Get the top SNP that matches PLINK formatting
    top_snp <- loclist[[cell]]$index_snp 

    print(top_snp)

    chr_num <- chr
    
    print(chr_num)
    
    # Construct the PLINK bash command
    plink_cmd <- paste0(
      # requires PLINK 1.9 (1.9.0-b.7.11) on PATH
      "plink --bfile ", genotype_files, "/chr", chr_num, " ",
      "--chr ", chr_num, " ",
      "--r2 --ld-snp ", top_snp, " ",
      "--ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 ",
      "--out ", plink_out, "/", target_gene, "_mat_eqtl_", cell
    )
    print(plink_cmd)

    message("Running PLINK for: ", cell)
    system(plink_cmd)
    
    ld_file <- paste0(plink_out, "/", target_gene, "_mat_eqtl_", cell, ".ld")
    
    if (file.exists(ld_file)) {
      plink_ld <- read_table(ld_file, show_col_types = FALSE)
      
      ld_to_merge <- plink_ld %>%
        dplyr::select(MarkerID = SNP_B, ld = R2)
      
      # Append the LD data to the locus object
      loclist[[cell]]$data <- loclist[[cell]]$data %>%
        dplyr::select(-any_of("ld")) %>% 
        dplyr::left_join(ld_to_merge, by = "MarkerID") 
      
      # not super sure if critical, but doing it anyway:
      #loclist[[cell]]$data$pos <- as.numeric(loclist[[cell]]$data$pos)
      message("Successfully merged LD for: ", cell)
    } else {
      warning("PLINK failed to create LD file for ", cell)
    }
  }

print("Appended LD, now save loclist as RDS object for plotting")

saveRDS(loclist, paste0(plink_out, "/tenk10k_", target_gene, "_locus_plot_data.RDS")) # this data can be plotted in R on local computer 