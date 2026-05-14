# Concatenate SMR results from a study

library(tidyverse)
library(arrow)
library(glue)
library(fs)
library(qvalue)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params
STUDY <- snakemake@wildcards[['study']]
# PARAMS <- list(
#     fdr_msmr = 0.05
# ) 

read_data <- function(x){
    file_parts <- path_split(x)[[1]]
    biosample <- file_parts[4]
    phenotype <- file_parts[5]

    data <- read_tsv_arrow(x) %>%
        mutate(biosample = biosample, phenotype = phenotype) |> 
        relocate(biosample, phenotype)
    return(data)
}

calc_q <- function(p, ...) {
    possibly(qvalue, qvalue_truncp(p))(p, ...) %>% .$qvalues
}

INPUT <- list(
    msmr = dir_ls(glue("results/smr/{STUDY}"), recurse = TRUE, glob = "*.msmr"),
    snps = dir_ls(glue("results/smr/{STUDY}"), recurse = TRUE, glob = "*.snps4msmr.list")
)

df_msmr <- map_df(INPUT[['msmr']], read_data)  %>% 
    mutate(qval_msmr_biosample_pheno = calc_q(p_SMR_multi), .by = c("biosample", "phenotype")) %>% 
    mutate(qval_msmr_biosample = calc_q(p_SMR_multi), .by = c("biosample")) %>% 
    mutate(qval_msmr_pheno = calc_q(p_SMR_multi), .by = c("phenotype")) %>% 
    mutate(qval_msmr = calc_q(p_SMR_multi))

# write results
write_parquet(df_msmr, OUTPUT[['msmr']], compression = "gzip")

# read and concatenate snps 4 msmr results
df_snps <- map_df(INPUT[['snps']], read_data)

write_parquet(df_snps, OUTPUT[['snps']], compression = "gzip")


# qvalue adjustment

# write results

# df_all <- df_all %>% 
#     group_by(biosample) %>% 
#     mutate(qval_msmr = qvalue(p_SMR_multi)$qvalues,
#            qval_singlesmr = qvalue(p_SMR)$qvalues)