#!/usr/bin/env bash
# Build the public bundle of summary-level figure inputs for deposit on Zenodo.
# Run via `make bundle` from the repository root after `make link-data`, which
# writes INVENTORY.tsv first; symbolic links are followed, so the archive holds
# the files themselves:
#
#   bash make_figure_bundle.sh <repo>/config/figure_inputs.txt out.tar.gz [--include-review]
#
# Only files explicitly listed in the manifest are packed (no directory globs).
# Lines tagged REVIEW (cell-level data) are skipped unless --include-review is
# given, after confirming they may be shared. Raw single-cell data, genotypes and
# read files are always refused. The archive root holds README.md (from
# BUNDLE_README.md), INVENTORY.tsv and figure_inputs.sha256; the build stops if
# the inventory's checksums disagree with the files.
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

[[ -f INVENTORY.tsv ]] || { echo "INVENTORY.tsv not found; run make inventory" >&2; exit 1; }
readme="scripts/0-preprocess/BUNDLE_README.md"

shasum -a 256 "${files[@]}" > figure_inputs.sha256
# the inventory must list exactly the bundled files, with the same checksums
if ! diff <(awk '{print $2 "\t" $1}' figure_inputs.sha256 | sort) \
          <(awk -F'\t' -v inc="$include_review" 'NR > 1 && (inc == "--include-review" || $9 != "yes") {print $1 "\t" $3}' INVENTORY.tsv | sort) >&2; then
    echo "INVENTORY.tsv does not match the files being bundled; rerun make inventory" >&2
    rm figure_inputs.sha256
    exit 1
fi
if [[ "$include_review" != "--include-review" ]]; then
    awk -F'\t' 'NR == 1 || $9 != "yes"' INVENTORY.tsv > INVENTORY.bundle.tsv
else
    cp INVENTORY.tsv INVENTORY.bundle.tsv
fi
tar -czhf "$out" --transform "s|^${readme}\$|README.md|; s|^INVENTORY.bundle.tsv\$|INVENTORY.tsv|" \
    "$readme" INVENTORY.bundle.tsv figure_inputs.sha256 "${files[@]}"
rm figure_inputs.sha256 INVENTORY.bundle.tsv
echo "wrote $out with ${#files[@]} file(s); $skipped REVIEW file(s) skipped"
