
library(data.table)
library(tidyverse)

pheno <- snakemake@wildcards[["pheno"]]
gwas <- snakemake@input[["gwas"]]
trait_metadata <- snakemake@input[["trait_metadata"]]
output <- snakemake@output[[1]]
threads  <- snakemake@threads

# interactive mode
# pheno <- "crohns"
# gwas <- "resources/sumstats/gwas/crohns.gwas"
# trait_metadata <- "resources/metadata/trait_metadata_curated.xlsx"
# output <- "resources/pipeline_ma/crohns.ma"
# threads <- 8
setDTthreads(threads)

# Summary stats from https://doi.org/10.1038/s41588-023-01384-0
# ndex variant chosen as the most significant variant in the locus and annotated as CHR:POS:A1:A2. 
# CHR, chromosome; POS, genomic position in genome build 38; A1, reference allele; A2, effect allele. 
# bOR and P-value are from the inverse-variance-weighted fixed-effect meta-analysis (two-tailed) 
# including all EAS samples. cNearest gene to the index variant. EA, effect allele; EAF, effect allele frequency.

df_trait_meta <- readxl::read_xlsx(trait_metadata)
setDT(df_trait_meta)

# Make alleles uppercase
# Add N value
# Select those that are only in non-finnish europeans
# No need to liftover, we are using hg38

gwas_df <- fread(gwas) %>% 
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2),
           CHR = as.numeric(CHR)) %>%
    filter(!is.na(AF_NFE), CHR %between% c(1,22))


# match with tenk10k variants
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

gwas_df[df_tenk10k_freq,
        `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2),
        on = c("CHR" = "CHR", "BP" = "POS")]

# Remove variants with missing alleles (unmatched CHR:BP)
gwas_df <- gwas_df[!is.na(A1) & !is.na(A2)]

# match allele based on alphabetical order
# suffix: .x for gwas_df, .y for tenk10k
gwas_df[, `:=`(Amin.x = pmin(Allele1, Allele2),
               Amax.x = pmax(Allele1, Allele2),
               Amin.y = pmin(A1, A2),
               Amax.y = pmax(A1, A2))]

# get estimated number of effective sample size
n_eff <- df_trait_meta[trait_id == pheno, n_eff]

df_out <- gwas_df %>% 
    # filter variant with allele mismatch
    filter(Amin.x == Amin.y & Amax.x == Amax.y) %>% 
    mutate(b = fifelse(Allele2 == A1, BETA_NFE, -BETA_NFE),
            # this is using the original AF_NFE, assuming Allele2 is faithfully reflecting AF_NFE
           freq = fifelse(Allele2 == A1, AF_NFE, 1-AF_NFE),
           n = n_eff) %>% 
    # filter based on maf >1%
    filter(pmin(freq, 1 - freq) > 0.01) %>%
    select(SNP = snp_id, A1, A2, freq, b,
           se = SE_NFE, p = P_NFE, n)

fwrite(df_out, output, sep = "\t")
