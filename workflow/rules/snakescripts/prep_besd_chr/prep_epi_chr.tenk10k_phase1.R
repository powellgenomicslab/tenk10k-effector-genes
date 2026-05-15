## Purpose: Prepare .epi files for BESD input from GTF annotation and cell-specific data
## Input: GTF file, BESD source directory
## Output: .epi files per cell and chromosome
## NOTE: Hardcoded paths below (/g/data/fy54/, /g/data/ei56/) are specific to NCI Gadi; update for other environments
library(rtracklayer)
library(data.table)
library(fs)
library(tidyverse)
library(glue)

GTF <- "/g/data/fy54/reference/GRCh38-gencode-v44/genes/genes.gtf.gz"
SOURCE_DIR <- "/g/data/ei56/as8574/analysis/TenK10K_SMR/inputs/besd"

CELLS <- dir_ls(SOURCE_DIR) %>% basename()

df_gtf <- readGFF(GTF)
setDT(df_gtf)
df_gtf[, seqid := gsub("^chr", "", seqid)]

for (c in CELLS) {
    for (chr in 1:22) {
        EPI <- glue("{SOURCE_DIR}/{c}/{c}_Chr{chr}.epi")
        df_epi <- fread(EPI)
        df_epi[df_gtf[type == "gene"], `:=`(start = i.start, end = i.end),
               on = c("V2" = "gene_id")]
        df_epi[, mid := as.integer(start + (floor((end - start) / 2)))]

        df_out <- df_epi[, .(V1, V2, V3, mid, V5, V6)]

        OUTDIR <- glue("resources/besd/tenk10k_phase1/{c}")
        OUTFILE <- glue("{OUTDIR}/chr{chr}.epi")
        dir_create(OUTDIR)
        fwrite(df_out, OUTFILE, sep="\t", quote=F, col.names=F)
    }
}