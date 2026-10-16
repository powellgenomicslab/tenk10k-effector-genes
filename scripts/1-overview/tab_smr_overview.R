source("scripts/0-preprocess/preprocess_results.R")
source("scripts/util/write_table.R")

# average hits per disease / cell type / combination, broken down by MR sensitivity criteria.
#
# `sig` doesn't exist under preprocess_strict.R (only the older preprocess.R defined it) --
# every use below is `mr` instead, matching this pipeline's evidence-flag columns
# (mr, sensitivity, coloc, mr_sens, mr_coloc, mr_sens_coloc; see
# workflow/rules/snakescripts/preprocess/tenk10k_phase1.v5.R for how these are derived,
# including the NA-propagation fix in commit 378d8da).
df_summary <- df_msmr_tenk10k %>%
  mutate(pheno_label = fct_reorder(pheno_label, mr, sum),
         trait_type = ifelse(supercategory == "biological", "Biomarker", "Disease")) %>%
  group_by(pheno_label, cell_type, trait_type, pheno_cat) %>%
  summarise(
    n_mr = sum(mr, na.rm = TRUE),
    n_single_instrument = sum(n_snps_smr == 1, na.rm = FALSE),
    median_n_instrument_msmr = median(n_snps_smr[n_snps_smr > 1], na.rm = FALSE),
    n_heidi_pass = sum(mr & p_HEIDI >= 0.05, na.rm = TRUE),
    n_het_ivw_pass = sum(mr & phet_ivw >= 0.05, na.rm = TRUE),
    n_pleio_mrlink2_pass = sum(mr & psigmay_mrlink2 >= 0.05, na.rm = TRUE),
    n_sens = sum(mr_sens, na.rm = TRUE),
    n_coloc_pph4_pass = sum(mr & coloc_pph4 >= 0.8, na.rm = TRUE),
    n_mvcoloc_pph4_pass = sum(mr & mvcoloc_pph4 >= 0.8, na.rm = TRUE),
    n_coloc = sum(mr_coloc, na.rm = TRUE),
    n_sens_coloc = sum(mr_sens_coloc, na.rm = TRUE),
    n_test = n(),
    .groups = "drop"
  ) %>%
  mutate(
    prop_mr = n_mr / n_test,
    prop_sens = n_sens / n_test,
    prop_coloc = n_coloc / n_test,
    prop_sens_coloc = n_sens_coloc / n_test
  ) %>%
  arrange(desc(pheno_label), desc(n_mr))

write_table(df_summary, "mr_summary", 3)
