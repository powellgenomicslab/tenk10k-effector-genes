# Purpose: Prepare Locus zoom of eQTL

# Interactive job
# qsub -I -q normal -P fy54 -l mem=32GB,storage=gdata/ei56+gdata/io72+gdata/fy54,ncpus=1
# conda deactivate
# project=smr/smr_tenk10k
# cd /g/data/ei56/rt3501/${project}
# module use /g/data/io72/apps/Modules/modulefiles
# module load R/4.4.3 # example R

####################### Read me for locuszoomr  #######################

# https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html 
# The genomic locus can be specified in several ways. The simplest is to specify a gene by name/symbol using the gene argument. The location of the gene is obtained from the specified Ensembl database. 
# The amount of flanking regions can either be controlled by specifying flank which defaults to 50kb either side of the ends of the gene. 
# flank can either be a single number or a vector of 2 numbers if different down/upstream flanking lengths are required. 
# Alternatively a fixed genomic window (eg. 1 Mb) centred on the gene of interest can be specified using the argument fix_window. The locus can be specified manually by specifying the chromosome using seqname and genomic position range using xrange. 
# Finally, a region can be specified by naming the index_snp, in which case the object data is searched for the coordinates of that SNP and the size of the region defined using fix_window or flank.
# You need an access token emailed to you to use the LDlinkR API (token = "3593e031216b"). However, they require SNPs in either chrCHR:BP format or rsID format, and they use the 1000 Genomes reference which may not have your SNPs.
# Therefore, it is easier to simply run PLINK using the genotype files on your index SNPs

# By default, running the locus() function assigns a SNP as the index SNP (access it through locus$index_snp)
# Calculate the r2 for all SNPs in the region to the index SNP
# You can run PLINK within the Rscript (rather than doing it separately in a bash script)

# You can explicitly label the SNP used for MR and not indicate the index SNP. 

####################### Load packages and connect db for locuszoomr functions #######################
setwd("/g/data/fy54/rt3501/repos/tenk10k-causal")

library(locuszoomr)
library(tidyverse)
library(data.table)
library(ragg)
library(cowplot)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(LDlinkR)

# Fetch your annotation database, ensDb_v106 is more updated then v86 but most of the protein coding genes are the same. 
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
gene <- "GPX1"
gene_ens <- "ENSG00000233276"
cellType <- c("CD4_TCM", "cDC2", "ILC")
topeqtlsnps <- mr %>% dplyr::filter(Gene == gene & biosample %in% cellType) %>% dplyr::select(biosample, topSNP) %>% deframe() %>% as.list()
chr <- mr %>% dplyr::filter(Gene == gene) %>% dplyr::pull(ProbeChr) %>% unique()

topeqtlsnps
gene
gene_ens 
gwas_name
cellType
####################### Read in GWAS data for locuszoomr #######################
gwas <- fread(gwas_fp) 
# need to add the chr, bp columns back for locus_plot(). fyi this is suuuuper slow so maybe use a datatable version?
gwas <- gwas %>%
  tidyr::separate(SNP, into = c("chrom", "pos"), sep = ":", extra = "drop", remove = FALSE) %>%
  dplyr::rename(MarkerID = SNP) # So it matches the eQTL results 

print("Read in GWAS file & make sure SNP is now MarkerID")
####################### Set the gene symbol and ensembl ID, chromsome and (later) cell type of interest #######################

# Create list of the eQTL files across all cell types, you need to filter out your GOI
#eqtl_fp <- paste0(matrix_eqtl_folders, "/", cellType, "/", "common_raw.tsv")
eqtl_fp <- paste0(eqtl_folders, "/", cellType, "/", "chr", chr, "_meqtl.tsv")

locus_input <- lapply(eqtl_fp, fread)

# for SAIGE eQTL
# locus_input <- locus_input %>%
#   map(~ .x %>% dplyr::rename(p = "p.value") %>%
#   dplyr::filter(gene == gene_ens))

# for "Matrix eQTL" formatted files, input: 
locus_input <- locus_input %>%
  map(~ .x %>% tidyr::separate(SNP, into = c("chrom", "pos"), sep = ":", extra = "drop", remove = FALSE) %>%
  dplyr::rename(p = "p-value", MarkerID = SNP) %>%
  dplyr::filter(gene == gene_ens))

names(locus_input) <- cellType

print("Read in eQTL files")

locus_input[["GWAS"]] <- gwas

# Run locus on everything in locus_input (NAMED LIST)
loclist <- lapply(locus_input, locus, gene = gene, ens_db = "EnsDb.Hsapiens.v86", flank = 5e5)

# save(gwas_loc, loclist, gene, max_y, file = "resources/crohns_case_study/figures/locus_zoom_plot_objects.RData")

############### GET PLINK DATA 

# Loop through every item in your locus plot data list
for (cell in names(loclist)) {
    # Get the top SNP that matches PLINK formatting
    top_snp <- loclist[[cell]]$index_snp 
    
    print(top_snp)
    # Extract the chromosome number from the SNP 
    chr_num <- stringr::str_split(top_snp, ":")[[1]][1]
    print(chr_num)
    
    # Construct the PLINK bash command
    plink_cmd <- paste0(
      "module use /g/data/io72/apps/Modules/modulefiles && ",
      "module load plink/1.9.0-b.7.11 && ",
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
saveRDS(loclist, paste0(figures, gene, "_loclist.RDS"))
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
# print(mr_stats)

# final_plot_names = c("Crohns Disease", 
#   "CD4 TCM (P MR = 4.37e-7, pHEIDI = 0.01, \COLOC PP H4 = 0.96)",
#   "cDC2 (P MR = 1.59e-11, pHEIDI = 0.67, \nCOLOC PP H4 = 0.96)", 
#   "ILC (Not tested; not an eGene)"
# )
final_plot_names = c("Crohns Disease", 
  "CD4 TCM - MR Associated",
  "cDC2 - MR Associated", 
  "ILC - Not tested; not an eGene"
)
png(paste0(figs, gene, "_GPX1_locus_plot_CD4TCM_cDC2.png"), width = 6, height = 15, units = "in", res = 300, type = "cairo") # Note it's res, not dpi. 

#png(paste0(figs, gene, "_multi-celltype_locus_zoom_with_LD_matrixeqtl.png"), width = 6, height = 15, units = "in", res = 300, type = "cairo") # Note it's res, not dpi. 
max_y = 30
multi_layout(nrow = 4,
             plots = {
               locus_plot(loclist[[4]], use_layout = FALSE, legend_pos = 'topleft', pcutoff = 5e-08, labels = c(topeqtlsnps[["CD4_TCM"]], topeqtlsnps[["cDC2"]]), highlight = gene, main = final_plot_names[1], ylim = c(0, max_y))
               locus_plot(loclist[[1]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topeqtlsnps[["CD4_TCM"]], highlight = gene, main = final_plot_names[2], ylim = c(0, max_y))
               locus_plot(loclist[[2]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topeqtlsnps[["cDC2"]], highlight = gene, main = final_plot_names[3], ylim = c(0, max_y))
               locus_plot(loclist[[3]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = NULL, highlight = gene, main = final_plot_names[4], ylim = c(0, max_y))
             })
dev.off()

svg(paste0(figs, gene, "_multi-celltype_locus_zoom_with_LD_matrixeqtl.svg"), width = 6, height = 15) # Note it's res, not dpi. 
max_y = 30
multi_layout(nrow = 4,
             plots = {
               locus_plot(loclist[[4]], use_layout = FALSE, legend_pos = 'topleft', pcutoff = 5e-08, labels = c(topeqtlsnps[["CD4_TCM"]], topeqtlsnps[["cDC2"]]), highlight = gene, main = final_plot_names[1], ylim = c(0, max_y))
               locus_plot(loclist[[1]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topeqtlsnps[["CD4_TCM"]], highlight = gene, main = final_plot_names[2], ylim = c(0, max_y))
               locus_plot(loclist[[2]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = topeqtlsnps[["cDC2"]], highlight = gene, main = final_plot_names[3], ylim = c(0, max_y))
               locus_plot(loclist[[3]], use_layout = FALSE, legend_pos = NULL, pcutoff = 5e-08, labels = NULL, highlight = gene, main = final_plot_names[4], ylim = c(0, max_y))
             })
dev.off()











# Try ggplot - This works, BUT LESS PRETTY 
# max_y = 20
# library(cowplot)
# p1 <- locus_ggplot(loclist[[4]], legend_pos = 'topleft', pcutoff = 5e-08, labels = my_labels, ylim = c(0, max_y), highlight = gene) + ggtitle(final_plot_names[1])
# p2 <- locus_ggplot(loclist[[1]], legend_pos = NULL, pcutoff = 5e-08, labels = my_labels, ylim = c(0, max_y), highlight = gene) + ggtitle(final_plot_names[2])
# p3 <- locus_ggplot(loclist[[2]], legend_pos = NULL, pcutoff = 5e-08, labels = my_labels, ylim = c(0, max_y), highlight = gene) + ggtitle(final_plot_names[3])
# p4 <- locus_ggplot(loclist[[3]], legend_pos = NULL, pcutoff = 5e-08, labels = my_labels, ylim = c(0, max_y), highlight = gene) + ggtitle(final_plot_names[4])

# p_all <- plot_grid(p1, p2, p3, p4, nrow = 4)
# ggsave(paste0("resources/crohns_case_study/figures/", gene, "_multi-celltype_gg_locus_with_LD_matrixeqtl.png"), p_all, device = ragg::agg_png(),
#        width = 6, height = 15, units = "in", bg = "white", scaling = 1, dpi = 300)

# save(loclist, file = paste0("resources/crohns_case_study/figures/", gene, "_multi-celltype_locus_zoom_plot_objects_matrixeqtl.RData"))
