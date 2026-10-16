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

# Nimbus Sans (URW base 35 fonts, AGPL-3.0 with a font exception): a free font
# with the metrics of Helvetica, which the figure scripts ask for (as Helvetica or
# Arial) but cannot be redistributed. conda-forge does not package it, and without
# it fontconfig substitutes the wider DejaVu Sans, so text overflows the panels.
# Install the pinned release with URW's fontconfig rule mapping Helvetica to Nimbus
# Sans, plus the same mapping for Arial and Arial Narrow. Typst does not read
# fontconfig, so the Typst rules pass --font-path $CONDA_PREFIX/fonts.
urw=20200910
curl -fsSL -o "$tmp/urw.tar.gz" "https://github.com/ArtifexSoftware/urw-base35-fonts/archive/refs/tags/${urw}.tar.gz"
tar -xzf "$tmp/urw.tar.gz" -C "$tmp"
mkdir -p "$CONDA_PREFIX/fonts" "$CONDA_PREFIX/etc/fonts/conf.d"
cp "$tmp/urw-base35-fonts-${urw}"/fonts/NimbusSans*.otf "$CONDA_PREFIX/fonts/"
cp "$tmp/urw-base35-fonts-${urw}/fontconfig/urw-nimbus-sans.conf" "$CONDA_PREFIX/etc/fonts/conf.d/61-urw-nimbus-sans.conf"
cat > "$CONDA_PREFIX/etc/fonts/conf.d/62-arial-nimbus-sans.conf" <<'CONF'
<?xml version="1.0"?>
<!DOCTYPE fontconfig SYSTEM "fonts.dtd">
<!-- written by figures.post-deploy.sh: Arial has Helvetica's metrics -->
<fontconfig>
  <alias binding="same"><family>Arial</family><accept><family>Nimbus Sans</family></accept></alias>
  <alias binding="same"><family>Arial Narrow</family><accept><family>Nimbus Sans Narrow</family></accept></alias>
</fontconfig>
CONF
fc-cache -f "$CONDA_PREFIX/fonts"
