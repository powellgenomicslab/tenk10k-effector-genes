library(data.table)
library(dplyr)
library(readxl)

setDTthreads(snakemake@threads)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PHENO <- snakemake@wildcards[['phenotype']]

n_eff <- read_excel(INPUT$trait_metadata) %>% 
    filter(include) %>%
    filter(trait_id == PHENO) %>% 
    pull(n_eff)
df_ma <- fread(INPUT$ma) %>% 
    mutate(n = n_eff) %>% 
    select(Predictor = SNP, A1, A2, n,
           Direction = b, P = p) %>% 
    unique() %>% 
    # filter all missing Direction (LDAK doesn't allow missing information)
    filter(!is.na(Direction), !is.na(P)) %>% 
    filter(!duplicated(Predictor) & !duplicated(Predictor, fromLast = TRUE),
           A1 != A2)

fwrite(df_ma, OUTPUT[[1]], sep = "\t", na = "NA", quote = FALSE)