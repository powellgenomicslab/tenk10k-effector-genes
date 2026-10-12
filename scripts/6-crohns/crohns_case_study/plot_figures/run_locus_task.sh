#!/bin/bash

# Arguments
GENE=$1
Chr=$2
PROBE=$3

# # Directories
# dir=/g/data/ei56/rt3501
# # go to where you relative paths in Rscript are set to
# cd ${dir}
# logs=${dir}/scripts/figures/1-generate_locus_zoom_object_${GENE}.log
# # conda activate tidyverse -> don't need this, using modules refering to renv with packages pre-installed
# Rscript smr/smr_onek1k/smr_analysis_CD_Liu2023/scripts/figures/1-generate_locus_zoom_object.R ${GENE} ${Chr} 2>&1 >> ${logs}

PROJECT_DIR="/g/data/fy54/rt3501/repos/tenk10k-causal"
SCRIPT_DIR="${PROJECT_DIR}/scripts/crohns_case_study/plot_figures"
LOG_DIR="${PROJECT_DIR}/scripts/crohns_case_study/logs/figures"  # Defining a dedicated log prefix

RSCRIPT_FILE="${SCRIPT_DIR}/1-generate_locus_zoom_object.R"
LOG_FILE="${LOG_DIR}/1-generate_locus_zoom_object_${GENE}.log"

mkdir -p "${LOG_DIR}"

cd "${PROJECT_DIR}"
Rscript "${RSCRIPT_FILE}" "${GENE}" "${Chr}" "${PROBE}">> "${LOG_FILE}" 2>&1

