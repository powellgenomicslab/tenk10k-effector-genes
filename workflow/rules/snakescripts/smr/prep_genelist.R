library(rtracklayer)
library(data.table)
library(fs)
library(tidyverse)

GTF <- snakemake@input[[1]]
OUTPUT <- snakemake@output[[1]]

# GTF <- "resources/smr_misc/tenk10k_phase1.gtf.gz"
# OUTPUT <- "resources/smr_misc/tenk10k_phase1.genelist.txt"

df_gtf <- readGFF(GTF)

setDT(df_gtf)

df_gtf[, chr := gsub("^chr", "", seqid) %>% as.numeric()]

df_gtf_subset <- df_gtf %>% 
    filter(type == "gene" & str_detect(seqid, "^chr")) %>%
    mutate(chr = gsub("^chr", "", seqid) %>% as.numeric()) %>% 
    filter(!is.na(chr)) %>%
    arrange(chr, start, end) %>% 
    select(chr, start, end, gene_id, strand)

fwrite(df_gtf_subset, OUTPUT, sep="\t", quote=F, col.names=F)
