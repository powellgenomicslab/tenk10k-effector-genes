#!/usr/bin/env bash
# Build the public bundle of summary-level figure inputs for deposit on
# Zenodo / Hugging Face. Run on Gadi from the data root (the tenk10k-causal
# working directory), pointing at this checkout's manifest:
#
#   bash make_figure_bundle.sh <repo>/config/figure_inputs.txt out.tar.gz [--include-review]
#
# Only files explicitly listed in the manifest are packed (no directory globs).
# Lines tagged REVIEW (cell-level data, in-sample LD, third-party files) are
# skipped unless --include-review is given, after confirming they may be shared.
# Raw single-cell data, genotypes and read files are always refused.
set -euo pipefail

manifest="${1:?usage: make_figure_bundle.sh <figure_inputs.txt> <out.tar.gz> [--include-review]}"
out="${2:?usage: make_figure_bundle.sh <figure_inputs.txt> <out.tar.gz> [--include-review]}"
include_review="${3:-}"

raw_pattern='\.(h5ad|h5|loom|mtx|bam|cram|fastq|fq|bed|bim|fam|pgen|pvar|psam|bgen|vcf)(\.gz)?$|(^|/)genotypes/|(^|/)matrix_eqtl/'

files=()
skipped=0
while IFS= read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    path="${line%%$'\t'*}"
    # tracked repository files are not bundled
    [[ "$path" == metadata/* || "$path" == scripts/* || "$path" == config/* ]] && continue
    if [[ "$path" =~ $raw_pattern ]]; then
        echo "refusing raw-data path: $path" >&2
        exit 1
    fi
    if [[ "$line" == *"# REVIEW:"* && "$include_review" != "--include-review" ]]; then
        echo "skipping (REVIEW): $path" >&2
        skipped=$((skipped + 1))
        continue
    fi
    [[ -f "$path" ]] || { echo "missing: $path" >&2; exit 1; }
    files+=("$path")
done < "$manifest"

shasum -a 256 "${files[@]}" > figure_inputs.sha256
tar -czhf "$out" figure_inputs.sha256 "${files[@]}"
rm figure_inputs.sha256
echo "wrote $out with ${#files[@]} file(s); $skipped REVIEW file(s) skipped"
