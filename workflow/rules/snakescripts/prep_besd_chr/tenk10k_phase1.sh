#!/bin/bash -l
# Purpose: Prepare BESD and EPI files for TenK10K phase1 analysis
# Input: BESD source directory, GTF file
# Output: Linked BESD/ESI files and processed EPI files per cell and chromosome
#
# This script stages pre-built BESD/ESI files; it does not build them. The
# SAIGE-QTL -> BESD conversion lives in the external TenK10K_SMR repository
# (see workflow/README.md). Override the two paths below via the environment,
# e.g. TENK10K_BESD_DIR=/my/besd GENCODE_GTF=/my/genes.gtf.gz.
set -euo pipefail

# prepare BESD file for tenk10k phase1

SOURCE_DIR="${TENK10K_BESD_DIR:?set TENK10K_BESD_DIR to the directory of pre-built BESD files}"
GTF="${GENCODE_GTF:?set GENCODE_GTF to the GENCODE v44 gene annotation GTF}"

CELLS=($(find "${SOURCE_DIR}" -mindepth 1 -maxdepth 1 -type d | xargs basename -a))

for CELL in "${CELLS[@]}"; do
    echo "Processing cell: ${CELL}"
    
    # Create output directory for the cell
    OUTDIR="resources/besd/tenk10k_phase1/${CELL}"
    mkdir -p "${OUTDIR}"

    # Prepare BESD file for the cell
    for CHR in {1..22}; do
        # link BESD and ESI files for each chromosome
        ln -s "${SOURCE_DIR}/${CELL}/${CELL}_Chr${CHR}.besd" "${OUTDIR}/chr${CHR}.besd"
        ln -s "${SOURCE_DIR}/${CELL}/${CELL}_Chr${CHR}.esi" "${OUTDIR}/chr${CHR}.esi"

    done
done

# Process EPI
conda activate renv
Rscript --vanilla "workflow/rules/snakescripts/prep_besd_chr/prep_epi_chr.tenk10k_phase1.R"
