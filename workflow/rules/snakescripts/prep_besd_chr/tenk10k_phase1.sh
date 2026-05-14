#!/bin/bash -l

# prepare BESD file for tenk10k phase1

SOURCE_DIR="/g/data/ei56/as8574/analysis/TenK10K_SMR/inputs/besd"
GTF="/g/data/fy54/reference/GRCh38-gencode-v44/genes/genes.gtf.gz"

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
