# scripts to preprocess GWAS summary statistics for FnBmd
# data: https://qlu-lab.org/data.html
# publication: https://www.nature.com/articles/s41588-024-01934-0

# Columns in original file:
# CHR = chromosome
# BP = base pair position (hg19 genome build)
# SNP = rsID 
# A1 = effect allele (i.e., count of this allele is used in GWAS)
# A2 = non-effect allele
# EAF = A1 allele frequency
# BETA = effect size
# SE = standard error
# Z = z-score
# P = p-value
# N_eff = effective sample size

library(data.table)
library(tidyverse)

setDTthreads(snakemake@threads)

liftover_script <- snakemake@input[["liftover_script"]]
gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# gwas_file <- "resources/sumstats/gwas/FnBmd.gwas"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"
# output_filepath <- "resources/pipeline_ma/FnBmd.ma"

source(liftover_script)
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

# Read in summary statistics
gwas_df <- fread(gwas_file)[CHR %in% 1:22]

# Liftover code
coord_df <- gwas_df[, .(CHR, BP, SNP)]
coord_df <- coord_df[, CHR := paste0("chr", CHR)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("CHR", "pos_b38", "SNP")

out_df <- merge(gwas_df, coord_hg38_df, by = c("CHR", "SNP"))

# match with tenk10k variants

out_df[df_tenk10k_freq,
        `:=`(snp_id = i.SNP, A1.y = i.A1, A2.y = i.A2),
        on = c("CHR" = "CHR", "pos_b38" = "POS")]

out_df <- out_df[!is.na(snp_id)]

out_df[, `:=`(Amin.x = pmin(A1, A2),
               Amax.x = pmax(A1, A2),
               Amin.y = pmin(A1.y, A2.y),
               Amax.y = pmax(A1.y, A2.y))]


out_df <- out_df %>% 
    # filter variant with allele mismatch
    filter(Amin.x == Amin.y & Amax.x == Amax.y) %>% 
    mutate(b = fifelse(A1 == A1.y, BETA, -BETA),
           freq = fifelse(A1 == A1.y, EAF, 1-EAF)) %>% 
    # filter based on maf >1%
    filter(pmin(freq, 1 - freq) > 0.01) %>%
    select(SNP = snp_id, A1 = A1.y, A2 = A2.y, freq, b,
           se = SE, p = P, n = N)

output_filepath <- snakemake@output[[1]]
fwrite(out_df, output_filepath, sep = "\t")
