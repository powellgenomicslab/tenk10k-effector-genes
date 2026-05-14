
library(arrow)

OUTPUT <- snakemake@output
df_magma <- read_parquet("results/aggregate/tenk10k_phase1.magma.parquet.gz")
df_msmr <- read_parquet("results/aggregate/tenk10k_phase1.msmr.parquet.gz")

# filter and recalculate results based on available genes in both MAGMA and TenK10K MSMR
gene_universe <- intersect(df_magma$GENE, df_msmr$probeID)

writeLines(gene_universe, OUTPUT$gene_universe)
