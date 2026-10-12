#!/usr/bin/env bash
# Run by Snakemake (inside the activated environment) after it creates the
# figures environment.
set -euo pipefail
: "${CONDA_PREFIX:?run inside the activated figures environment}"
rlib="$CONDA_PREFIX/lib/R/library"

# Use only this environment's packages: a personal R library (e.g.
# ~/Library/R/<arch>/4.4/library) holds binaries built for another R
# installation, which crash conda R when loaded.
mkdir -p "$CONDA_PREFIX/lib/R/etc"
printf '%s\n' "# written by figures.post-deploy.sh: ignore personal R libraries" \
    'R_LIBS_USER=""' >> "$CONDA_PREFIX/lib/R/etc/Renviron.site"

# sciscales (scientific-notation labels, used by scripts/7-drug/otp_combined.R)
# is only on GitHub; install the commit the published figures were made with.
# Fetched with the system curl and installed with R CMD INSTALL, avoiding
# remotes and the R curl package.
sha=8d68bab84d21fe0e36f271358f0f0eddbe87dbcd
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
curl -fsSL -o "$tmp/sciscales.tar.gz" "https://github.com/christyray/sciscales/archive/${sha}.tar.gz"
R CMD INSTALL --library="$rlib" "$tmp/sciscales.tar.gz"

# oaqc (orbit-aware quad census; used by the "backbone" layout in graphlayouts,
# scripts/5-pathways/enrichment_features.R) is on CRAN but not conda. Pin the
# version the figures were made with; fall back to the CRAN archive once a newer
# release replaces it.
oaqc=oaqc_2.0.0.tar.gz
curl -fsSL -o "$tmp/$oaqc" "https://cran.r-project.org/src/contrib/$oaqc" ||
    curl -fsSL -o "$tmp/$oaqc" "https://cran.r-project.org/src/contrib/Archive/oaqc/$oaqc"
R CMD INSTALL --library="$rlib" "$tmp/$oaqc"
