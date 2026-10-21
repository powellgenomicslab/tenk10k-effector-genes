#!/usr/bin/env bash
# Download the public bundle of summary-level inputs for the manuscript figures
# and tables, verify every file against its checksum, and unpack it into the
# repository root (results/ and resources/). Run via `make data`.
#
# If the record's files are restricted, a reader who has been granted access
# can download them with a Zenodo personal access token: set ZENODO_TOKEN, or
# save the token in ~/.zenodo_token (or the file named by ZENODO_TOKEN_FILE).
# The download then goes through Zenodo's API; without a token it is public.
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

token="${ZENODO_TOKEN:-}"
token_file="${ZENODO_TOKEN_FILE:-$HOME/.zenodo_token}"
if [[ -z "$token" && -r "$token_file" ]]; then
    token=$(tr -d '[:space:]' < "$token_file")
fi

# https://zenodo.org/records/<id>/files/<name>?download=1 -> the same file through the API
re='^https://(sandbox\.)?zenodo\.org/records/([0-9]+)/files/([^?]+)'
if [[ -n "$token" && "$url" =~ $re ]]; then
    api="https://${BASH_REMATCH[1]}zenodo.org/api/records/${BASH_REMATCH[2]}/files/${BASH_REMATCH[3]}/content"
    echo "Downloading with your Zenodo token: $api"
    # the header is read from stdin, so the token does not appear in the process list
    printf 'Authorization: Bearer %s\n' "$token" | curl -fL --retry 3 -H @- -o "$tmp/bundle.tar.gz" "$api"
else
    curl -fL --retry 3 -o "$tmp/bundle.tar.gz" "$url"
fi

# restricted files are answered with Zenodo's web page, not the archive
if ! gzip -t "$tmp/bundle.tar.gz" 2> /dev/null; then
    if head -c 512 "$tmp/bundle.tar.gz" | grep -qi '<html'; then
        echo "Zenodo returned a web page instead of the archive: the record's files are restricted." >&2
        if [[ -n "$token" ]]; then
            echo "Your token was used; check that your Zenodo account has been granted access to the record." >&2
        else
            echo "If you have been granted access, create a personal access token at zenodo.org" >&2
            echo "(Account > Applications) and save it in ~/.zenodo_token, or set ZENODO_TOKEN." >&2
        fi
    else
        echo "The download is not a valid gzip archive." >&2
    fi
    exit 1
fi
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
