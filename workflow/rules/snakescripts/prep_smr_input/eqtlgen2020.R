# get probe list and p-threshold for SMR from eQTLGen 2020

library(data.table)
library(fs)
library(readr)

setDTthreads(4)
input_file <- "resources/download/eqtl/eqtlgen_vosa2020/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
df <- fread(input_file)

p_thresh <- df[, max(Pvalue)]
probe_list <- df[, sort(unique(Gene))]

# write to output
outdir <- snakemake@output[[1]]
dir_create(path(outdir, "bulk_wb"))

write_lines(p_thresh, path(outdir, "bulk_wb", "pthresh.txt"))
write_lines(probe_list, path(outdir, "bulk_wb", "probe.txt"))