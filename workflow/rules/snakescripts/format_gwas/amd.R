# process AMD GWAS summary statistics
# source: Alex Hewitt

library(data.table)
library(tidyverse)
library(dbplyr)
library(DBI)
library(RSQLite)

pheno <- snakemake@wildcards[["pheno"]]
gwas <- snakemake@input[["gwas"]]
trait_metadata <- snakemake@input[["trait_metadata"]]
output <- snakemake@output[[1]]
threads  <- snakemake@threads

# interactive mode
# pheno <- "amd"
# gwas <- "resources/sumstats/gwas/amd.gwas"
# trait_metadata <- "resources/metadata/trait_metadata_curated.xlsx"
# ref_db <- "resources/ensembl_vcf/ensembl_r115_merged.db"
# output <- "resources/pipeline_ma/amd.ma"
# threads <- 8
# setDTthreads(threads)

df_trait_meta <- readxl::read_xlsx(trait_metadata)
setDT(df_trait_meta)

# Make alleles uppercase
# Add N value
# Select those that are only in non-finnish europeans
# No need to liftover, we are using hg38

df_gwas <- fread(gwas, key = "MarkerName") %>% 
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))

# get chr pos from SNP rsID based on sqlite ensembl vcf

# 1. Connect to your new database
con <- dbConnect(SQLite(), ref_db)
db_var <- tbl(con, "variants")

# 1. Drop the old simple index (it's redundant now)
# dbExecute(con, "DROP INDEX IF EXISTS idx_rsid")

# 2. Create the Covering Index
# We include 'chrom' and 'pos' in the index so SQLite never has to look at the main table
# dbExecute(con, "CREATE INDEX idx_rsid_cover ON variants(rsid, chrom, pos)")

# Optimize connection for speed
dbExecute(con, "PRAGMA synchronous = OFF")
dbExecute(con, "PRAGMA journal_mode = MEMORY")
dbExecute(con, "PRAGMA cache_size = 1000000") # Uses ~1GB RAM for cache

# 2. Load your input data (8M rows)
total_rows <- nrow(df_gwas)
chunk_size <- 50000 # 10k rows is a sweet spot for SQLite

# Create indices for splitting
df_gwas[, split_id := ceiling(seq_len(.N) / chunk_size)]
n_chunks <- max(df_gwas$split_id)

query_db <- function(i, idx = "MarkerName") {
  df_chunk <- df_gwas[split_id == i]
  copy_to(con, df_chunk, "temp_chunk", indexes = list(idx), temporary = TRUE, overwrite = TRUE)
  db_chunk <- tbl(con, "temp_chunk")

  db_chunk %>%
    left_join(tbl(con, "variants"), by = c("MarkerName" = "rsid")) %>%
    select(everything(), chrom, pos) %>%
    collect()
}

list_results <- list()

for(i in 1:n_chunks) {

  # D. Store result
  list_results[[i]] <- query_db(i)

  # Progress bar
  cat(sprintf("\rProcessed chunk %d / %d (%.1f%%)", i, n_chunks, (i/n_chunks)*100))
}

df_annot <- rbindlist(list_results) |> 
  # remove variants with missing chr/pos
  filter(!is.na(chrom) & !is.na(pos)) |> 
  mutate(chrom = as.integer(chrom),
         pos = as.integer(pos))

# match with tenk10k variants
df_tenk10k_freq <- fread("resources/genotypes_frq/tenk10k_phase1.frq")

df_annot[df_tenk10k_freq,
        `:=`(snp_id = i.SNP, A1 = i.A1, A2 = i.A2),
         on = c("chrom" = "CHR", "pos" = "POS")]

# match allele based on alphabetical order
# suffix: .x for gwas_df, .y for tenk10k
df_annot[, `:=`(Amin.x = pmin(Allele1, Allele2),
               Amax.x = pmax(Allele1, Allele2),
               Amin.y = pmin(A1, A2),
               Amax.y = pmax(A1, A2))]

# get estimated number of effective sample size
calc_n_eff <- function(n_case, n_control) {
    4 / (1/n_case + 1/n_control)
}

# from accompanying metadata
n_eff <- calc_n_eff(
    n_case = 105345,
    n_control = 1232600
)

df_out <- df_annot %>% 
    # filter variant with allele mismatch
    filter(Amin.x == Amin.y & Amax.x == Amax.y) %>% 
    mutate(b = fifelse(Allele1 == A1, Effect, -Effect),
            # this is using the original AF_NFE, assuming Allele2 is faithfully reflecting AF_NFE
           freq = fifelse(Allele1 == A1, Freq1, 1-Freq1),
           n = n_eff) %>% 
    # filter based on maf >1%
    # filter(pmin(freq, 1 - freq) > 0.01) %>%
    select(SNP = snp_id, A1, A2, freq, b,
           se = StdErr, p = `P-value`, n)

fwrite(df_out, output, sep = "\t")
