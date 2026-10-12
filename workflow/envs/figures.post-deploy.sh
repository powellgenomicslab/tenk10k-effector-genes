#!/usr/bin/env bash
# Run by Snakemake after it creates the figures environment.
# sciscales (scientific-notation labels, used by scripts/7-drug/otp_combined.R)
# is only on GitHub; install the commit the published figures were made with.
# The tarball is fetched with the system curl and installed with R CMD INSTALL,
# avoiding remotes/R curl, which can crash in conda R on macOS.
set -euo pipefail
sha=8d68bab84d21fe0e36f271358f0f0eddbe87dbcd
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
curl -fsSL -o "$tmp/sciscales.tar.gz" "https://github.com/christyray/sciscales/archive/${sha}.tar.gz"
R CMD INSTALL "$tmp/sciscales.tar.gz"
