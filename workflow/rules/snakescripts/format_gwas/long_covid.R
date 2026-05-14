# format long covid data
# source: https://www.nature.com/articles/s41588-025-02100-w
# sumstats from: https://my.locuszoom.org/gwas/192226/
# use Strict case definition vs. broad control
# alt = effect allele

library(data.table)
library(tidyverse)

N_case <- 3018
N_control <- 994582
setDTthreads(snakemake@threads)

# Read in summary statistics
gwas_df <- fread(snakemake@input[[1]]) %>% 
    rename(CHR = "#chrom", POS = "pos") %>% 
    filter(CHR %in% 1:22) %>% 
    mutate(CHR = as.numeric(CHR))

# merge with TenK10K frequency data
tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

# align alleles and select columns to output
output_df <- gwas_df %>%
    merge(tenk10k_freq, by = c("CHR", "POS")) %>% 
    filter(!is.na(SNP)) %>% 
    mutate(beta = fifelse(alt == A1, beta, -beta),
           freq = fifelse(alt == A1, alt_allele_freq, 1 - alt_allele_freq),
           N = floor(4 / (1/N_case + 1/N_control)),
           p = 10^(-neg_log_pvalue)) %>% 
    select(SNP, A1, A2, freq, b = beta, se = stderr_beta, p, N)

# Write results    
fwrite(output_df, snakemake@output[[1]], sep = "\t", na = "NA", quote = FALSE)