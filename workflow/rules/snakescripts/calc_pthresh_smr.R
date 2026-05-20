# calculate p-threshold for SMR per celltype

library(data.table)

INPUT <- snakemake@input[[1]]
OUTPUT <- snakemake@output[[1]]

df <- fread(INPUT)

df_thresh <- df[, .(pthresh = max(top_pval)), by = celltype]

fwrite(df_thresh, OUTPUT, sep = "\t")