# Processing summary statistics from Grove et al. 2019 (ASD)
# DOI: https://doi.org/10.1038/s41588-019-0344-8
#
# TenK10K (tenk10k_phase1.frq) is treated as the reference panel for alleles
# and allele frequencies — the output .ma file is harmonised onto the TenK10K
# A1 so that beta and freq are consistent with the eQTL BESD/ESI used by SMR.
#
# Input  (iPSYCH-PGC_ASD_Nov2017): CHR SNP BP A1 A2 INFO OR SE P
#   A1 = effect allele (OR is reported w.r.t. A1); A2 = other allele; build hg19.
#   No allele frequency is provided in the source sumstats.
# Output (.ma for SMR):            SNP A1 A2 freq b se p n
#   A1 = effect allele = TenK10K reference A1; A2 = other allele; build hg38.
#   freq = TenK10K allele frequency of A1 (from tenk10k_phase1.frq), NOT a
#   source-GWAS frequency (none available).
#
# Harmonisation steps:
#   1. Liftover SNP coordinates hg19 -> hg38.
#   2. Join to TenK10K .frq by (CHR, pos_b38) to pick up reference alleles
#      (A1_ref, A2_ref) and the TenK10K freq of A1_ref (A1_freq).
#   3. Drop variants whose unordered allele pair disagrees with the reference
#      (catches multi-allelic / strand mismatches).
#   4. Re-orient the GWAS effect (log OR) onto A1_ref: flip sign when the GWAS
#      effect allele (A1) is the reference's other allele. The output `freq`
#      is taken directly from A1_freq (already w.r.t. A1_ref); no flip needed.

library(data.table)
library(tidyverse)

setDTthreads(snakemake@threads)

pheno             <- snakemake@wildcards[["pheno"]]
gwas_file         <- snakemake@input[["gwas"]]
chain_file        <- snakemake@input[["hg19tohg38"]]
liftover_script   <- snakemake@input[["liftover_script"]]
trait_metadata    <- snakemake@input[["trait_metadata"]]

# interactive mode
# pheno           <- "asd"
# gwas_file       <- "resources/sumstats/gwas/asd.gwas"
# chain_file      <- "resources/misc/hg19ToHg38.over.chain"
# liftover_script <- "workflow/rules/snakescripts/hg19tohg38.R"
# trait_metadata  <- "resources/metadata/trait_metadata_curated.xlsx"

df_trait_meta <- readxl::read_xlsx(trait_metadata)
setDT(df_trait_meta)

gwas_df <- fread(gwas_file) %>%
    filter(CHR %in% 1:22)
gwas_df[, CHR := as.numeric(CHR)]

# Liftover hg19 -> hg38
source(liftover_script)

coord_df <- gwas_df[, .(CHR, BP, SNP)]
coord_df[, CHR := paste0("chr", CHR)]
colnames(coord_df) <- c("seqnames", "start", "snpid")

coord_hg38_df <- liftover2hg38(coord_df, chain_file)
colnames(coord_hg38_df) <- c("CHR", "pos_b38", "SNP")

out_df <- merge(gwas_df, coord_hg38_df, by = c("CHR", "SNP"))

# Annotate with TenK10K alleles and allele frequency
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

out_df[df_tenk10k_freq,
       `:=`(snp_id = i.SNP, A1_ref = i.A1, A2_ref = i.A2, A1_freq = i.MAF),
       on = c("CHR", "pos_b38" = "POS")]

n_eff <- df_trait_meta[trait_id == pheno, n_eff]

out_df_format <- out_df %>%
    filter(!is.na(snp_id)) %>%
    # Alphabetical allele sort to detect mismatches between GWAS and reference
    mutate(Amin_gwas = pmin(A1, A2),
           Amax_gwas = pmax(A1, A2),
           Amin_ref  = pmin(A1_ref, A2_ref),
           Amax_ref  = pmax(A1_ref, A2_ref)) %>%
    filter(Amin_gwas == Amin_ref & Amax_gwas == Amax_ref) %>%
    # Orient effect to TenK10K A1; A1_freq is already freq(A1_ref) from the .frq join
    mutate(b    = ifelse(A1 == A1_ref, log(OR), -log(OR)),
           freq = A1_freq,
           n    = n_eff) %>%
    filter(pmin(freq, 1 - freq) > 0.01) %>%
    select(SNP = snp_id, A1 = A1_ref, A2 = A2_ref, freq, b,
           se = SE, p = P, n)

out_file <- snakemake@output[[1]]
# out_file <- "resources/pipeline_ma/asd.ma"
fwrite(out_df_format, out_file, sep = "\t", na = "NA", quote = FALSE)
