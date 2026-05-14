# update magma output with hgnc symbol and fdr

library(tidyverse)
library(data.table)
library(qvalue)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

# Read in magma genes.out file
df <- fread(INPUT$out)
df_gene <- fread(INPUT$gene_loc, header = FALSE) %>% 
    select(gene_id = V1, hgnc_symbol = V6)

df[df_gene, `:=`(gene = i.hgnc_symbol),
    on = c("GENE" = "gene_id")]

df[, `:=`(qval = qvalue(P)$qvalues,
          lfdr = qvalue(P)$lfdr,
          p_bh = p.adjust(P, method = "BH"))]

fwrite(df, OUTPUT[[1]], sep = "\t")