#!/usr/bin/env bash
# Download the public bundle of summary-level inputs for the manuscript figures
# and tables, verify every file against its checksum, and unpack it into the
# repository root (results/ and resources/). Run via `make data`.
#
# usage: download_figure_bundle.sh <bundle_url>
set -euo pipefail

url="${1:?usage: download_figure_bundle.sh <bundle_url>}"
if [[ "$url" == *TODO* ]]; then
    echo "The figure-input bundle has not been deposited yet (bundle_url in config/figures.yaml)." >&2
    exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

curl -fL --retry 3 -o "$tmp/bundle.tar.gz" "$url"
mkdir "$tmp/unpacked"
# unpack into an empty scratch directory first; refuse absolute or parent paths
if tar -tzf "$tmp/bundle.tar.gz" | grep -qE '^/|(^|/)\.\.(/|$)'; then
    echo "Refusing bundle with absolute or parent-directory paths." >&2
    exit 1
fi
tar -xzf "$tmp/bundle.tar.gz" -C "$tmp/unpacked"

if command -v shasum > /dev/null; then sha256check="shasum -a 256 -c"; else sha256check="sha256sum -c"; fi
(cd "$tmp/unpacked" && $sha256check figure_inputs.sha256)

(cd "$tmp/unpacked" && find results resources -type f 2>/dev/null || true) | while read -r f; do
    mkdir -p "$(dirname "$f")"
    cp "$tmp/unpacked/$f" "$f"
done
echo "Figure inputs installed under results/ and resources/."
