args <- commandArgs(trailingOnly = TRUE)
#args <- c("lymph", "B_memory")

trait <- args[1]
celltype <- args[2]
output_filename <- args[3]
smr_filepath <- args[4]
snp_smr_filepath <- args[5]

library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)
library(genetics.binaRies)
library(arrow)
library(tibble)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(qvalue)


setDTthreads(7)

# Define paths
## Mount /g/data1a/ei56/as8574/analysis/TenK10K_SMR to /workspace if running on NCI
## Mount /g/data/fy54/analysis/tenk10k-causal to /data if running on NCI
data_dir <- "/data"
genotype_dir <- "/genotypes"
work_dir <- "/workspace"

# Point to PLINK inside singularity
plink_path <- "/usr/local/bin/plink"

check_harmonisation <- function(harmonised_dat, bim_df){
  allele_check <- bim_df[SNP %in% toupper(harmonised_dat$SNP)]
  harmonised_allele_order <- as.data.table(harmonised_dat[, c("SNP", 
  "effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", 
  "other_allele.outcome")])
  harmonised_allele_order$SNP <- toupper(harmonised_allele_order$SNP)
  allele_check <- merge(harmonised_allele_order, allele_check, by = "SNP")
  eaf_exposure_check <- allele_check$effect_allele == allele_check$`effect_allele.exposure`
  eaf_outcome_check <- allele_check$effect_allele == allele_check$`effect_allele.outcome`
  
  return(all(eaf_exposure_check) & all(eaf_outcome_check))
}

prep_exposure <- function(x, y){
  exposure <- y[["gene"]]
  exp_dat <- format_data(as.data.frame(x), type = "exposure",
    phenotype_col = "phenotype",
    snp_col = "variant", 
    gene_col = "gene",
    beta_col = "beta",
    se_col = "se",
    chr_col = "chr",
    pos_col = "pos",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "p-value"
  )
  return(exp_dat)
}

run_ivw_mr <- function(x, y, outcome_dat = NULL, ld_dir = NULL){
    exp_dat <- prep_exposure(x, y)
    harmonised_dat <- harmonise_data(exp_dat, outcome_dat, action = 2)

    # Extract SNPs, read in LD matrix
    snp_list <- harmonised_dat$SNP
    
    # For error handling
     empty_row <- data.frame(
          n_snps = length(snp_list),
          estimate = NA,
          se = NA,
          lower_ci = NA,
          upper_ci = NA,
          pval = NA,
          residual_se = NA,
          het_tstat = NA,
          het_pval = NA,
          fstat = NA,
          i_sq = NA,
          max_ld = NA,
          note = NA
      )

    # Get LD matrix
    # Gene-specific LD files
    # PLINK format: Chr,Chr:POS:REF:ALT,Dist,POS,ALT,REF
    bim_filepath <- file.path(ld_dir, paste0(y[["gene"]], ".bim"))
    ld_filepath <- file.path(ld_dir, paste0(y[["gene"]], ".ld"))

    if (! file.exists(bim_filepath) | ! file.exists(ld_filepath)){
        error_message <- paste0("No LD matrix or BIM for ", y[["gene"]])
        message(error_message)
        return(empty_row %>% mutate(note = error_message))
    }

    bim_df <- fread(bim_filepath, col.names = c("ProbeChr", "SNP", "dist", "POS", "other_allele", "effect_allele"))
    ld_mat <- as.matrix(as.data.frame(fread(ld_filepath, header = FALSE)))
    colnames(ld_mat) <- tolower(bim_df$SNP)
    rownames(ld_mat) <- tolower(bim_df$SNP)
    gene_ld_mat <- ld_mat[snp_list, snp_list]

    # Empty row
    if (!is.matrix(gene_ld_mat)){
      # Run IVW MR with correlation matrix
      mr_obj <- dat_to_MRInput(harmonised_dat)
          
      ## Append LD matrix to the MRInput object
      mr <- mr_ivw(mr_obj[[1]], model = "default")
      max_correlation <- NA

      # Extract values
      result <- data.frame(
          n_snps = mr@SNPs,
          estimate = mr@Estimate,
          se = mr@StdError,
          lower_ci = mr@CILower,
          upper_ci = mr@CIUpper,
          pval = mr@Pvalue,
          residual_se = mr@RSE,
          het_tstat = mr@Heter.Stat[1],
          het_pval = mr@Heter.Stat[2],
          fstat = mr@Fstat,
          i_sq = (mr@Heter.Stat[1] - (mr@SNPs - 1))/mr@Heter.Stat[1],
          max_ld = max_correlation,
          note = NA

      )
      return(result)
    }

    if (check_harmonisation(harmonised_dat, bim_df)){
      # Run IVW MR with correlation matrix
      mr_obj <- dat_to_MRInput(harmonised_dat)
          
      ## Append LD matrix to the MRInput object
      mr_obj[[1]]@correlation <- gene_ld_mat
      mr <- mr_ivw(mr_obj[[1]], correl = TRUE, model = "random")

      # Calulate highest LD value
      tested_ld_matrix <- mr_obj[[1]]@correlation
      max_correlation <- max((tested_ld_matrix[upper.tri(tested_ld_matrix)]^2))

      # Extract values
      result <- data.frame(
          n_snps = mr@SNPs,
          estimate = mr@Estimate,
          se = mr@StdError,
          lower_ci = mr@CILower,
          upper_ci = mr@CIUpper,
          pval = mr@Pvalue,
          residual_se = mr@RSE,
          het_tstat = mr@Heter.Stat[1],
          het_pval = mr@Heter.Stat[2],
          fstat = mr@Fstat,
          i_sq = (mr@Heter.Stat[1] - (mr@SNPs - 1))/mr@Heter.Stat[1],
          max_ld = max_correlation,
          note = NA

      )
      return(result)
    } else{
        error_message <- paste0("Alleles are incorrectly harmonised for ", y[["gene"]])
        message(error_message)
        return(empty_row %>% mutate(note = error_message))
    }
}

test_chrom <- function(chrom = NULL, msmr_tenk10k = NULL, snp_smr_df = NULL){
  print(chrom)
  # Chromosome-specific files
  ld_dir <- file.path(data_dir, sprintf("resources/ld/eqtl/tenk10k_phase1/chr%d", chrom))
  exposure_path <- sprintf("/data/resources/matrix_eqtl/tenk10k_phase1/%s/chr%d_meqtl.tsv", celltype, chrom)
  exposure_snp_path <- sprintf("/data/resources/besd/tenk10k_phase1/%s/chr%d.esi", celltype, chrom)
  outcome_path <- sprintf("/data/resources/ma_by_chr/%s/chr%d.ma", trait, chrom)
  bfile_filepath <- sprintf("/genotypes/chr%d_common_variants", chrom)

  # Prepare outcome data
  outcome_dat <- read_outcome_data(
    filename = outcome_path,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "N"
  )

  # Narrow down SMR info
  smr_df <- msmr_tenk10k[ProbeChr == chrom]
  snp_df <- snp_smr_df[gene %in% smr_df$probeID]

  # Merge SNP-Gene pairs
  # Read in eQTL data
  eqtl_df <- fread(exposure_path, sep = "\t", header = TRUE)
  esi_df <- fread(exposure_snp_path, sep = "\t", header = FALSE, col.names = c("chr", "SNP", "dist", "pos", "effect_allele", "other_allele", "eaf"))

  # Extract data for SNP-Gene pairs in eQTL data and calculate SE
  combined_df <- merge(snp_df, eqtl_df, by.x = c("gene", "variant"), by.y = c("gene", "SNP"), all.x = TRUE)
  combined_df <- merge(combined_df, esi_df, by.x = "variant", by.y = "SNP", all.x = TRUE)
  combined_df[, se := abs(beta)/abs(`t-stat`)]
  combined_df[, phenotype := trait]

  # Trim columns to only have useful information
  combined_df <- combined_df[, .(phenotype, variant, gene, chr, pos, effect_allele, other_allele, eaf, beta, se, `p-value`)]

  # Convert to tibble
  combined_df <- as_tibble(combined_df)

  # Run IVW MR for each gene
  result_df <- combined_df %>%
    group_by(gene) %>%
    group_modify(~ run_ivw_mr(.x, .y, outcome_dat, ld_dir = ld_dir)) %>%
    ungroup()

  # Combine with mSMR to compare
  output_df1 <- smr_df %>% select(biosample, phenotype, probeID, ProbeChr, 
    Gene, Probe_bp, p_eQTL, p_SMR, p_SMR_multi, p_HEIDI, nsnp_HEIDI)

  # Combine with IVW-LD results
  output_df <- output_df1 %>%
    left_join(result_df, by = c("probeID" = "gene"))

  return(output_df)
}

# Read in SMR file
msmr_tenk10k <- read_parquet(smr_filepath, as_data_frame = FALSE) |>
  filter(phenotype == trait, biosample == celltype, lfdr_msmr_pheno < 0.05) |>
  as.data.table()

# Read in snp list used in SMR
snp_smr_df <- read_parquet(snp_smr_filepath, as_data_frame = FALSE) |> 
  filter(biosample == celltype & phenotype == trait & gene %in% msmr_tenk10k$probeID) |>
  as.data.table()


chrom_list <- unique(msmr_tenk10k$ProbeChr)
ivw_chrom_list <- lapply(chrom_list, function(x) test_chrom(chrom = x, msmr_tenk10k =  msmr_tenk10k, snp_smr_df = snp_smr_df))

ivw_df <- bind_rows(ivw_chrom_list)
ivw_df <- ivw_df %>% arrange(ProbeChr, Probe_bp)

# Write out results
fwrite(ivw_df, output_filename, sep = "\t")