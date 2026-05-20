# calculate p-threshold for SMR per celltype
# generate probe list for SMR per celltype

library(data.table)
library(fs)
library(tidyverse)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

df <- fread(INPUT[[1]])

cells <- unique(df$celltype)

process <- function(c) {
    dir_create(path(OUTPUT, c))
    thresh <- df[celltype == c, max(top_pval)]
    probe_list <- df[celltype == c, unique(gene)]
    # fwrite(df_thresh, path(OUTPUT, c, "pthresh.txt"), sep = "\t")
    write_lines(thresh, path(OUTPUT, c, "pthresh.txt"))
    write_lines(probe_list, path(OUTPUT, c, "probe.txt"))
}

walk(cells, process)