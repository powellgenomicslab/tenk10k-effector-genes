"""
Map Ensembl gene IDs / HGNC symbols → STRING protein IDs.

Reads a gencode-style gene annotation file and queries the STRING
``get_string_ids`` API in batches to build a persistent lookup table.

Input:  gencode TSV (columns: chr, start, end, ensembl_gene_id, hgnc_symbol, gene_type)
Output: TSV with columns:
          ensembl_gene_id, hgnc_symbol, stringId, preferredName

Only protein-coding genes are mapped (STRING covers proteins).
The output can be joined in downstream prep scripts (e.g. mr_strict.R)
to replace Ensembl IDs with STRING IDs before submission.
"""

import io
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# ── Snakemake interface ──────────────────────────────────────────────
INPUT  = snakemake.input
OUTPUT = snakemake.output
PARAMS = snakemake.params
LOG    = snakemake.log

# ── Logging ──────────────────────────────────────────────────────────
log_file = str(LOG[0])
os.makedirs(os.path.dirname(log_file), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_file, mode="w"),
        logging.StreamHandler(sys.stderr),
    ],
)
log = logging.getLogger(__name__)

# ── Parameters ───────────────────────────────────────────────────────
species        = int(getattr(PARAMS, "species", 9606))
string_version = str(getattr(PARAMS, "string_version", "12.0"))
caller_id      = str(getattr(PARAMS, "caller_identity", "tenk10k_smr"))
batch_size     = int(getattr(PARAMS, "batch_size", 2000))
retries        = int(getattr(PARAMS, "retries", 3))

STRING_API = "https://string-db.org/api"

log.info("Start: map gencode → STRING IDs")
log.info(f"  gencode:        {INPUT.gencode}")
log.info(f"  output:         {OUTPUT.mapping}")
log.info(f"  species:        {species}")
log.info(f"  string_version: {string_version}")
log.info(f"  api_url:        {STRING_API}")
log.info(f"  batch_size:     {batch_size}")


# ── API helper ───────────────────────────────────────────────────────

def string_post_tsv(endpoint, data, max_retries=retries, timeout=300):
    """POST to STRING API, return TSV text."""
    url = f"{STRING_API}/tsv/{endpoint}"
    for attempt in range(max_retries):
        time.sleep(1)  # STRING asks ≥1 s between calls
        try:
            resp = requests.post(url, data=data, timeout=timeout)
            resp.raise_for_status()
            return resp.text
        except Exception as exc:
            if attempt < max_retries - 1:
                wait = 2 ** (attempt + 1)
                log.warning(f"  Retry {attempt+1}/{max_retries}: {exc} (wait {wait}s)")
                time.sleep(wait)
            else:
                raise


# ── Load gencode ─────────────────────────────────────────────────────
t0 = time.time()

gencode = pd.read_csv(str(INPUT.gencode), sep="\t")
log.info(f"Loaded {len(gencode)} genes from gencode")

# Filter to protein-coding only (STRING maps proteins)
pc = gencode[gencode["gene_type"] == "protein_coding"].copy()
log.info(f"Protein-coding genes: {len(pc)}")

# We'll query by hgnc_symbol (preferred by STRING for human) and
# keep the ensembl_gene_id for joining back.  Drop rows without
# a usable symbol.
pc = pc[pc["hgnc_symbol"].notna() & (pc["hgnc_symbol"] != "")].copy()
pc = pc.drop_duplicates(subset="hgnc_symbol")
log.info(f"Unique HGNC symbols to map: {len(pc)}")

symbols = pc["hgnc_symbol"].tolist()


# ── Batch mapping ────────────────────────────────────────────────────

all_mapped = []
n_batches = (len(symbols) + batch_size - 1) // batch_size

for i in range(0, len(symbols), batch_size):
    batch = symbols[i : i + batch_size]
    batch_num = i // batch_size + 1
    log.info(f"  Batch {batch_num}/{n_batches}: mapping {len(batch)} symbols …")

    try:
        txt = string_post_tsv("get_string_ids", {
            "identifiers":   "\r".join(batch),
            "species":       species,
            "limit":         1,
            "echo_query":    1,
            "caller_identity": caller_id,
        })
        if txt.strip():
            df = pd.read_csv(io.StringIO(txt), sep="\t")
            all_mapped.append(df)
            log.info(f"  Batch {batch_num}: mapped {len(df)} identifiers")
        else:
            log.warning(f"  Batch {batch_num}: empty response")
    except Exception as exc:
        log.error(f"  Batch {batch_num}: FAILED — {exc}")


# ── Assemble output ──────────────────────────────────────────────────

if all_mapped:
    mapped = pd.concat(all_mapped, ignore_index=True)
    log.info(f"Total mapped: {len(mapped)} rows")

    # Keep only the columns we need, rename for clarity
    keep = {}
    if "queryItem" in mapped.columns:
        keep["queryItem"] = "hgnc_symbol"
    if "stringId" in mapped.columns:
        keep["stringId"] = "stringId"
    if "preferredName" in mapped.columns:
        keep["preferredName"] = "preferredName"

    result = mapped.rename(columns=keep)[list(keep.values())].copy()

    # Deduplicate: keep first (best) match per symbol
    result = result.drop_duplicates(subset="hgnc_symbol")

    # Join back ensembl_gene_id
    ensembl_lookup = pc[["ensembl_gene_id", "hgnc_symbol"]].drop_duplicates(
        subset="hgnc_symbol"
    )
    result = result.merge(ensembl_lookup, on="hgnc_symbol", how="left")

    # Reorder columns
    cols = ["ensembl_gene_id", "hgnc_symbol", "stringId", "preferredName"]
    result = result[[c for c in cols if c in result.columns]]

    log.info(f"Final mapping: {len(result)} genes")
    log.info(f"  Mapped {result['stringId'].notna().sum()} / {len(pc)} "
             f"protein-coding genes")
else:
    log.warning("No genes mapped — writing empty file")
    result = pd.DataFrame(columns=[
        "ensembl_gene_id", "hgnc_symbol", "stringId", "preferredName"
    ])

# ── Write output ─────────────────────────────────────────────────────
out_path = str(OUTPUT.mapping)
os.makedirs(os.path.dirname(out_path), exist_ok=True)
result.to_csv(out_path, sep="\t", index=False)

elapsed = time.time() - t0
log.info(f"Done: {out_path} ({elapsed:.1f}s)")
