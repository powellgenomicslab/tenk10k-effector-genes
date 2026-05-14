#!/usr/bin/env bash

set -Eeuxo pipefail

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

PREFIX_BFILE_CHR="${snakemake_params[prefix_bfile_chr]}"
PREFIX_BESD_CHR="${snakemake_params[prefix_besd_chr]}"
PREFIX_OUT_CHR="${snakemake_params[prefix_out_chr]}"

cleanup() {
    rm -f "${PREFIX_OUT_CHR}"*
}

trap cleanup EXIT ERR

for i in {1..22}; do
    smr --bfile "${PREFIX_BFILE_CHR}${i}" \
        --gwas-summary ${snakemake_input[ma]} \
        --beqtl-summary "${PREFIX_BESD_CHR}${i}" \
        --peqtl-smr ${snakemake_params[p_eqtl_thresh]} \
        --extract-probe ${snakemake_input[probe]} \
        --maf ${snakemake_params[maf]} \
        --diff-freq-prop ${snakemake_params[diff_freq_prop]} \
        --smr-multi \
        --thread-num ${snakemake[threads]} \
        --out "${PREFIX_OUT_CHR}${i}" \
        2>&1 | tee "${snakemake_log[0]}" 
done

# Merge the results
awk 'NR == FNR || FNR > 1' "${PREFIX_OUT_CHR}"*.msmr > ${snakemake_output[smr]}
awk 'NR == FNR || FNR > 1' "${PREFIX_OUT_CHR}"*.snp_failed_freq_ck.list > ${snakemake_output[fail]}
cat "${PREFIX_OUT_CHR}"*.prbregion4msmr.list > ${snakemake_output[region]}

awk '
  BEGIN {print "gene\tvariant"}
  NR==1 {gene=$1; next} 
  $1 == "end" {get_next=1; next} 
  get_next {gene=$1; get_next=0; next} 
  {print gene "\t" $1}
' "${PREFIX_OUT_CHR}"*.snps4msmr.list > ${snakemake_output[instrument]}

