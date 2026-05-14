# quantify open targets target-disease evidence

# Using arrow / dplyr
library(tidyverse)
library(arrow)

INPUT <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params

# read trait metadata
# INPUT <- list(
#     otp_drug_dir = "resources/nci/otp_output/25.12/drug_molecule",
#     otp_drug_mechanism_dir = "resources/nci/otp_output/25.12/drug_mechanism_of_action",
#     otp_evidence_chembl_dir = "resources/nci/otp_output/25.12/evidence_chembl"
# )


ds_drug_mechanism <- open_dataset(INPUT$otp_drug_mechanism_dir)
df_drug_mechanism <- ds_drug_mechanism %>%
    select(chemblIds, actionType, mechanismOfAction, targetType, targets) %>% 
    collect() %>% 
    unnest(chemblIds) %>% 
    unnest(targets) %>% 
    rename(drugId = chemblIds, targetId = targets) %>%
    distinct()

# get drug annotations
ds_drug <- open_dataset(INPUT$otp_drug_dir)
df_drug <- ds_drug %>%
    select(id, drugType, name,
        #    yearOfFirstApproval, # removed in 26.03
          maximumClinicalStage
        ) %>% 
    collect() %>% 
    rename(drugId = id) %>%
    distinct()

ds_evidence <- open_dataset(INPUT$otp_evidence_chembl_dir)
df_evidence <- ds_evidence %>%
    # filter(datasourceId == "chembl") %>% 
    select(drugId, targetId, diseaseId, score) %>%
    collect() %>% 
    distinct()

fs::dir_create(PARAMS$output_dir)
write_tsv(df_drug_mechanism, OUTPUT[["drug_mechanism"]])
write_tsv(df_evidence, OUTPUT[["drug_evidence"]])
write_tsv(df_drug, OUTPUT[["drug_molecule"]])
