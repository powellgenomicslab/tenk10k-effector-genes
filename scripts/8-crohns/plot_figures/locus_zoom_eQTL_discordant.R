# Purpose: Prepare Locus zoom of eQTL
# Produces resources/crohns_case_study/figures/JAZF1_data.RDS (JAZF1 eQTL in CD4 Naive and cDC2, with LD),
# which plot_figures/JAZF1_scatterplot.R plots as Supplementary Figure 17.
# Requires individual-level genotypes (PLINK) to compute LD.

# Interactive job
# qsub -I -q normal -P fy54 -l mem=50GB,storage=gdata/ei56+gdata/io72+gdata/fy54,ncpus=1
# module load R/4.4.3 # example R

####################### Read me for locuszoomr  #######################

# https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html 
# The genomic locus can be specified in several ways. The simplest is to specify a gene by name/symbol using the gene argument. The location of the gene is obtained from the specified Ensembl database. 
# The amount of flanking regions can either be controlled by specifying flank which defaults to 50kb either side of the ends of the gene. 
# flank can either be a single number or a vector of 2 numbers if different down/upstream flanking lengths are required. 
# Alternatively a fixed genomic window (eg. 1 Mb) centred on the gene of interest can be specified using the argument fix_window. The locus can be specified manually by specifying the chromosome using seqname and genomic position range using xrange. 
# Finally, a region can be specified by naming the index_snp, in which case the object data is searched for the coordinates of that SNP and the size of the region defined using fix_window or flank.
# You need an access token to use the LDlinkR API. However, they require SNPs in either chrCHR:BP format or rsID format, and they use the 1000 Genomes reference which may not have your SNPs.
# Therefore, it is easier to simply run PLINK using the genotype files on your index SNPs

# By default, running the locus() function assigns a SNP as the index SNP (access it through locus$index_snp)
# Calculate the r2 for all SNPs in the region to the index SNP
# You can run PLINK within the Rscript (rather than doing it separately in a bash script)

# You can explicitly label the SNP used for MR and not indicate the index SNP. 

####################### Load packages and connect db for locuszoomr functions #######################
# Run from the repository root (holding resources/ and results/).

library(locuszoomr)
library(tidyverse)
library(data.table)
library(ragg)
library(cowplot)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(LDlinkR)

# Fetch your annotation database, ensDb_v106 is more updated then v86 but most of the protein coding genes are the same. So if you're 
# Important to specify localHub = TRUE
# Note: I ran these commands in the login node (with internet) ONCE, where I did not use the localHub arg
ah <- AnnotationHub(localHub = TRUE)
query(ah, c("EnsDb", "Homo sapiens"))
ensDb_v106 <- ah[["AH100643"]]

####################### Paths #######################

figs <- "resources/crohns_case_study/figures/"
gwas_name <- "crohns"
gwas_fp <- paste0("resources/pipeline_ma/", gwas_name, ".ma")
#eqtl_folders <- "resources/saige_eqtl/tenk10k_phase1"
eqtl_folders <- "resources/matrix_eqtl/tenk10k_phase1"

genotype_files="resources/genotypes/tenk10k_phase1"
plink_out = "resources/crohns_case_study/figures/ld/"
print("Set paths")

####################### Vars #######################
# set your gene of interest (or region)
# gene <- "ZFP36L1"
# gene_ens <- "ENSG00000185650"
# chr <- 14
# cellType <- c("MAIT", "CD4_TCM", "B_intermediate")

# get the top eQTL SNPs to label 
mr <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_raw.RDS")
gene <- "JAZF1"
gene_ens <- mr %>% dplyr::filter(Gene == gene) %>% pull(probeID) %>% unique()

cellType <- c("CD4_Naive", "cDC2")
topeqtlsnps <- mr %>% dplyr::filter(Gene == gene & biosample %in% cellType) %>% dplyr::select(biosample, topSNP) %>% deframe() %>% as.list()
chr <- mr %>% dplyr::filter(Gene == gene) %>% dplyr::pull(ProbeChr) %>% unique()

topeqtlsnps
gene
gene_ens 
gwas_name
cellType
####################### Read in GWAS data for locuszoomr #######################
# gwas <- fread(gwas_fp) 
# # need to add the chr, bp columns back for locus_plot(). fyi this is suuuuper slow so maybe use a datatable version?
# gwas <- gwas %>%
#   tidyr::separate(SNP, into = c("chrom", "pos"), sep = ":", extra = "drop", remove = FALSE) %>%
#   dplyr::rename(MarkerID = SNP) # So it matches the eQTL results 

# print("Read in GWAS file & make sure SNP is now MarkerID")
####################### Set the gene symbol and ensembl ID, chromsome and (later) cell type of interest #######################

# Create list of the eQTL files across all cell types, you need to filter out your GOI
#eqtl_fp <- paste0(matrix_eqtl_folders, "/", cellType, "/", "common_raw.tsv")
eqtl_fp <- paste0(eqtl_folders, "/", cellType, "/", "chr", chr, "_meqtl.tsv")

locus_input <- lapply(eqtl_fp, fread)

# for SAIGE eQTL
# locus_input <- locus_input %>%
#   map(~ .x %>% dplyr::rename(p = "p.value") %>%
#   dplyr::filter(gene == gene_ens))

# for Matrix eQTL: 
locus_input <- locus_input %>%
  map(~ .x %>% 
  tidyr::separate(SNP, into = c("chrom", "pos"), sep = ":", extra = "drop", remove = FALSE) %>%
  dplyr::rename(p = "p-value", MarkerID = SNP) %>%
  dplyr::filter(gene == gene_ens) %>%
  dplyr::mutate(p_transform = -log10(p)*sign(beta)))

names(locus_input) <- cellType

print("Read in eQTL files")

# locus_input[["GWAS"]] <- gwas

# Run locus on everything in locus_input (NAMED LIST)
# loclist <- lapply(locus_input, locus, gene = gene, ens_db = "EnsDb.Hsapiens.v86", flank = 5e5)
loclist <- lapply(locus_input, locus, gene = gene, ens_db = "EnsDb.Hsapiens.v86", flank = 5e5, yvar = "p_transform")

# save(gwas_loc, loclist, gene, max_y, file = "resources/crohns_case_study/figures/locus_zoom_plot_objects.RData")

############### GET PLINK DATA 
specific_top_snp <- "7:28135367:T:C"
# Loop through every item in your locus plot data list
for (cell in names(loclist)) {
    # Get the top SNP that matches PLINK formatting
    #top_snp <- loclist[[cell]]$index_snp 
    top_snp <- specific_top_snp

    print(top_snp)
    # Extract the chromosome number from the SNP 
    chr_num <- stringr::str_split(top_snp, ":")[[1]][1]
    print(chr_num)
    
    # Construct the PLINK bash command
    plink_cmd <- paste0(
      # requires PLINK 1.9 (1.9.0-b.7.11) on PATH
      "plink --bfile ", genotype_files, "/chr", chr_num, " ",
      "--chr ", chr_num, " ",
      "--r2 --ld-snp ", top_snp, " ",
      "--ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 ",
      "--out ", plink_out, gene, "_mat_eqtl_", cell
    )
    print(plink_cmd)

    message("Running PLINK for: ", cell)
    system(plink_cmd)
    
    ld_file <- paste0(plink_out, gene, "_mat_eqtl_", cell, ".ld")
    
    if (file.exists(ld_file)) {
      plink_ld <- read_table(ld_file, show_col_types = FALSE)
      
      ld_to_merge <- plink_ld %>%
        dplyr::select(MarkerID = SNP_B, ld = R2)
      
      # Append the LD data to the locus object
      loclist[[cell]]$data <- loclist[[cell]]$data %>%
        dplyr::select(-any_of("ld")) %>% 
        dplyr::left_join(ld_to_merge, by = "MarkerID") 
      
      # not super sure if critical, but doing it anyway:
      loclist[[cell]]$data$pos <- as.numeric(loclist[[cell]]$data$pos)

      message("Successfully merged LD for: ", cell)
    } else {
      warning("PLINK failed to create LD file for ", cell)
    }
  }

print("Appended LD, now plot")

######################## APPENDED LD, NOW PLOT ######### SAME SNP 

#final_plot_names = c("Crohns Disease", "MAIT", "CD4 TCM", "B intermediate")
# names(loclist) <- final_plot_names

# Apparently can't add 'topleft' in legend_pos for every sigle one ???? 
# my_labels = c("index", "14:68787474:C:T")
# topSNP = "14:68787474:C:T"

# png(paste0(figs, gene, "_multi-celltype_locus_zoom_with_LD_matrixeqtl.png"), width = 6, height = 15, units = "in", res = 300, type = "cairo") # Note it's res, not dpi. 
# max_y = 20
# multi_layout(nrow = 4,
#              plots = {
#                locus_plot(loclist[[4]], use_layout = FALSE, legend_pos = 'topleft', pcutoff = 5e-08, labels = my_labels, highlight = gene, main = final_plot_names[1], ylim = c(0, max_y))
#                locus_plot(loclist[[1]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topSNP, highlight = gene, main = final_plot_names[2], ylim = c(0, max_y))
#                locus_plot(loclist[[2]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topSNP, highlight = gene, main = final_plot_names[3], ylim = c(0, max_y))
#                locus_plot(loclist[[3]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topSNP, highlight = gene, main = final_plot_names[4], ylim = c(0, max_y))
#              })
# dev.off()

######################## APPENDED LD, NOW PLOT ######### DIFFERENT SNPS
# STATS 
# mr_stats <- mr %>% 
#   dplyr::filter(mr_sens_coloc) %>%
#   filter(Gene == "GPX1") %>% 
#   select(Gene, cell_type, biosample, p_SMR_multi, p_HEIDI, psigmay_mrlink2, coloc_pph4, coloc_pph4)

# final_plot_names = c("Crohns Disease", 
#   "CD4 Naive (pSMRmulti = 4.37e- 7, pHEIDI = 0.0133, \np-sigma-y = 0.552, coloc_pph4 = 0.964)",
#   "cDC2 (pSMRmulti = 1.59e-11, pHEIDI = 0.667, \np-sigma-y = 1.000, coloc_pph4 = 0.964)", 
#   "ILC (Not tested; not an eGene)"
# )
saveRDS(loclist, "resources/crohns_case_study/figures/JAZF1_data.RDS")
