#!/usr/bin/env bash
# Link the figure inputs listed in config/figure_inputs.txt from an existing
# data tree (e.g. the tenk10k-causal working directory on Gadi) into this
# checkout. Files are linked one by one, never whole directories, so anything
# the figure scripts write stays in this checkout and the data tree is left
# untouched. Run via `make link-data DATA_DIR=...`.
#
# usage: link_figure_inputs.sh <data_dir>
set -euo pipefail

data_dir="${1:?usage: link_figure_inputs.sh <data_dir>}"
data_dir="$(cd "$data_dir" && pwd)"
manifest="config/figure_inputs.txt"

linked=0
missing=0
while IFS=$'\t' read -r path _; do
    [[ -z "$path" || "$path" == \#* ]] && continue
    # tracked repository files are used as they are
    [[ "$path" == metadata/* || "$path" == scripts/* || "$path" == config/* ]] && continue
    [[ -e "$path" && ! -L "$path" ]] && continue
    if [[ -e "$data_dir/$path" ]]; then
        mkdir -p "$(dirname "$path")"
        ln -sfn "$data_dir/$path" "$path"
        linked=$((linked + 1))
    else
        echo "missing: $path" >&2
        missing=$((missing + 1))
    fi
done < "$manifest"

echo "linked $linked file(s) from $data_dir; $missing missing"
