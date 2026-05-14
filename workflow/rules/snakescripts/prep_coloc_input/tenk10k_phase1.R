library(data.table)
library(fs)
library(fst)
library(tidyverse)

setDTthreads(as.numeric(Sys.getenv("NCPUS")))
threads_fst(as.numeric(Sys.getenv("NCPUS")))
OUTPUT <- snakemake@output

df_gene_loc <- fread("resources/misc/gencode.v44.gene_type.tsv")
dir_eqtl <- path("resources/saige_eqtl/tenk10k_phase1/")
df_sig_eqtl <- fread("resources/brenner/tenk10k_phase1/common_eqtl.tsv")

cells <- unique(df_sig_eqtl$celltype)

process <- function(c) {
    dir_create(path(OUTPUT, c))
    
    egenes <- df_sig_eqtl[celltype == c, unique(gene)]

    df_eqtl <- fread(path(dir_eqtl, c, "common_raw.tsv")) %>% 
        filter(gene %in% egenes) %>% 
        mutate(varbeta = SE^2,
               MAF = pmin(AF_Allele2, 1-AF_Allele2)) %>% 
        select(gene, chr = CHR, beta = BETA, varbeta, position = POS, snp = MarkerID, N, MAF)

    for (chr_num in unique(df_eqtl$chr)) {
        dir_create(path(OUTPUT, c), "common_egenes")
        df_eqtl_chr <- df_eqtl[chr == chr_num] %>% 
            select(-chr)
        write_fst(df_eqtl_chr, path(OUTPUT, c, "common_egenes", paste0("chr", chr_num, ".fst")))
    }
}

walk(cells, process)