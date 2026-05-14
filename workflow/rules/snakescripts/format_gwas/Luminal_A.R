library(data.table)
library(tidyverse)

setDTthreads(8)
# setDTthreads(snakemake@threads)

liftover_script <- snakemake@input[["liftover_script"]]
gwas_file <- snakemake@input[["gwas"]]
chain_file <- snakemake@input[["chain_file"]]
pheno <- snakemake@wildcards[["pheno"]]
out_file <- snakemake@output[[1]]

# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# chain_file <- "resources/misc/hg19ToHg38.over.chain"
source(liftover_script)

# Get MAF information from TenK10K panel

df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

# pheno <- "Luminal_A"
# gwas_file <- paste0("resources/sumstats/gwas/", pheno, ".gwas")
# out_file <- paste0("resources/pipeline_ma/", pheno, ".ma")

gwas_df <- fread(gwas_file) %>% 
    filter(CHR %in% 1:22)

# Liftover code
coord_df <- gwas_df[, .(CHR, BP, var_name)]
coord_df <- coord_df[, CHR := paste0("chr", CHR)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("CHR", "pos_b38", "var_name")

out_df <- merge(gwas_df, coord_hg38_df, by = c("CHR", "var_name"))


out_df[df_tenk10k_freq,
       `:=`(snp_id = i.SNP, A1.ref = i.A1, A2.ref = i.A2, A1_freq = i.MAF),
       on = c("CHR", "pos_b38" = "POS")]

# format
out_df_format <- out_df %>% 
    filter(!is.na(snp_id)) %>%
    arrange(CHR, pos_b38) %>% 
    mutate(b = ifelse(A1 == A1.ref, Z * SE, -(Z*SE)),
           N = floor(N)) %>% 
    select(SNP = snp_id, A1 = A1.ref, A2 = A2.ref,
           freq = A1_freq, b, se = SE,  p = P, N)

# Write results
fwrite(out_df_format, out_file, sep = "\t", na = "NA", quote = FALSE)

