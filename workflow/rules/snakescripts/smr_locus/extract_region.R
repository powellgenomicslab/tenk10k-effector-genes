library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(glue)
library(tidyverse)
library(fs)
library(sys)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params
X <- snakemake@wildcards

# X <- list(study = "tenk10k_phase1",
# biosample = "Treg", pheno = "crohns",
# probe = "ENSG00000185670")
# PARAMS <- list(probe_flank_kb = 500)
# INPUT <- list(
#     saige = glue("resources/saige_eqtl/{X$study}/{X$biosample}/common_raw.tsv"),
#         gtf = glue("resources/smr_misc/{X$study}.gtf.gz"),
#         ma = glue("resources/ma/{X$pheno}.ma"),
#     dir_geno_chr = glue("resources/genotypes/{X$study}/"),
#     dir_besd_chr = glue("resources/besd/{X$study}/{X$biosample}/")
# )

# load data
gtf <- rtracklayer::import(INPUT$gtf)
# extract probe information
# Extract the region of interest using GRanges

df_gene <- as.data.table(gtf) %>% 
    filter(type == "gene", gene_id == X$probe) %>% 
    select(gene_id, gene_name, seqnames, start, end) %>% 
    mutate(chr = str_remove(seqnames, "chr") %>% as.numeric)
chr <- df_gene$chr
start <- df_gene$start - PARAMS$probe_flank_kb * 1000
end <- df_gene$end + PARAMS$probe_flank_kb * 1000
region <- GRanges(seqnames = df_gene$seqnames, 
                  ranges = IRanges(start = start, end = end))


df_ma <- fread(INPUT$ma) %>% 
    filter(str_detect(SNP, paste0("^", chr, ":")))
df_ma[, chr := chr]
df_ma[, pos := str_extract(SNP, "(?<=:)[0-9]+(?=:)") %>% as.numeric]

# query probe Info
query_eqtl_probe <- function(probe, chr = df_gene$chr) {
    out_prefix <- file_temp()
    besd_prefix <- fs::path(INPUT$dir_besd_chr, glue("chr{chr}"))

    args <- c("smr",
     "--beqtl-summary", besd_prefix,
     "--query", 1,
     "--probe", probe,
     "--out", out_prefix
    )
    exec_wait(args)
    results <- fread(paste0(out_prefix, ".txt"))
    file_delete(paste0(out_prefix, ".txt"))
    
    return(results)
}
df_eqtl <- query_eqtl_probe(X$probe)

# df_saige <- fread(INPUT$saige) 


# extract region
gtf_region <- subsetByOverlaps(gtf, region)

df_ma_region <- df_ma %>% 
    filter(chr == df_gene$chr &
           pos %between% c(start, end)) %>% 
    select(variant_id = SNP, chr, pos, ea = A1, oa = A2, eaf = freq, b, se, p)

df_eqtl_region <- df_eqtl %>% 
    select(variant_id = SNP, chr = Chr, pos = BP, ea = A1, oa = A2,
           eaf = Freq, b, se = SE, p)


export(gtf_region, OUTPUT$gtf, format = 'gtf')
fwrite(df_ma_region, OUTPUT$gwas, sep = "\t", quote = FALSE)
fwrite(df_eqtl_region, OUTPUT$eqtl, sep = "\t", quote = FALSE)

# extract LD for eqtl variants
top_eqtl <- df_eqtl_region %>% 
    .[variant_id %in% df_ma_region$variant_id] %>% 
    .[which.min(p), variant_id]

plink_ld <- function(topsnp, chr, start, end) {
    out_prefix <- file_temp()
    bfile_prefix <- fs::path(INPUT$dir_geno_chr, glue("chr{chr}"))
    args <- c(
            "plink",
                "--bfile", bfile_prefix,
                "--r2",
                "--chr", chr,
                "--from-bp", start,
                "--to-bp", end,
                # "--snps", snps,
                "--silent",
                "--ld-window-r2", 0,
                "--ld-window", 1e6,
                "--ld-snp", topsnp,
                "--out", out_prefix
        )
    exec_wait(args)
    results <- fread(paste0(out_prefix, ".ld")) %>%
        mutate(across(starts_with("SNP"), as.character))
    file_delete(paste0(out_prefix, ".log"))
    file_delete(paste0(out_prefix, ".ld"))
    return(results)
}

df_ld <- plink_ld(top_eqtl, df_gene$chr, start, end)

fwrite(df_ld, OUTPUT$ld, sep = "\t", quote = FALSE)