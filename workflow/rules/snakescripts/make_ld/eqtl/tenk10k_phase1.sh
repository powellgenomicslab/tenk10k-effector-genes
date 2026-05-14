#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

CHR=${snakemake_wildcards[chr]}
STUDY=${snakemake_wildcards[study]}
GENE_LOC=${snakemake_input[gene_loc]}

EQTL_FILE="resources/brenner/tenk10k_phase1/common_eqtl.tsv"

EGENES=($(awk -F, 'NR >1 {print $1}' $EQTL_FILE | sort -u))

for G in "${EGENES[@]}"; do
  echo "Processing eGene: $G"

  CHR=$(awk -v gene="$G" '$4 == gene {print $1; exit}' $GENE_LOC)
  START=$(awk -v gene="$G" '$4 == gene {print $2; exit}' $GENE_LOC)
  END=$(awk -v gene="$G" '$4 == gene {print $3; exit}' $GENE_LOC)

  echo "CHR: $CHR, START: $START, END: $END"
  
done
