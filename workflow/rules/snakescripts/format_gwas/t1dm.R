library(data.table)
library(tidyverse)
setDTthreads(4)

# Summary stats from https://doi.org/10.1038/s41586-021-03552-w
# NOTE: GRCh38 - use RSIDs to map back to GRCh37
gwas_df <- fread("TenK10K_SMR_nci/data/source_files/gwas/GCST90014023_buildGRCh38.tsv")
gwas_df <- gwas_df[chromosome %in% 1:22]
gwas_df[, chromosome := as.numeric(chromosome)]

# match with tenk10k variants
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

gwas_df[df_tenk10k_freq,
        `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2),
        on = c("chromosome" = "CHR", "base_pair_location" = "POS")]

gwas_df <- gwas_df[!is.na(A1) & !is.na(A2)]

gwas_df[, `:=`(Amin.x = pmin(effect_allele, other_allele),
               Amax.x = pmax(effect_allele, other_allele),
               Amin.y = pmin(A1, A2),
               Amax.y = pmax(A1, A2))]


df_out <- gwas_df %>% 
    # filter variant with allele mismatch
    filter(Amin.x == Amin.y & Amax.x == Amax.y) %>% 
    mutate(b = fifelse(effect_allele == A1, beta, -beta),
           freq = fifelse(effect_allele == A1, effect_allele_frequency, 1-effect_allele_frequency),
           n = sample_size) %>% 
    # filter based on maf >1%
    filter(pmin(freq, 1 - freq) > 0.01) %>%
    select(SNP = snp_id,  A1, A2, freq, b,
           se = standard_error, p = p_value, n)

output_filepath <- snakemake@output[[1]]
fwrite(df_out, output_filepath, sep = "\t")