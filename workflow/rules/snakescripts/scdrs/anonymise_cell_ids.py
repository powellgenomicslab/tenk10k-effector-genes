"""Replace the sequencing-pool part of scDRS cell IDs with anonymous codes.

Cell IDs are <10x barcode>_<pool>, e.g. GGCACGTTCTCATGGA_S0109 or
GCATTAGTCGCCATAA_S0188_re1. Each distinct pool label is mapped one-to-one to a
code (P001, P002, ...) in random order, and the barcode is kept, so IDs stay
unique and the same cell has the same ID in every file. The mapping is written
to POOL_MAP, which is not shared, and reused if it already exists so that
reruns give the same IDs. A pool missing from an existing map is an error.

Run by rule scdrs_anonymise_cell_ids in workflow/rules/scdrs.smk.
"""

import os

import numpy as np
import pandas as pd

INPUTS = list(snakemake.input.cells)
OUTPUTS = list(snakemake.output.cells)
POOL_MAP = snakemake.params.pool_map


def split_id(cell_id):
    barcode, _, pool = cell_id.partition("_")
    return barcode, pool


frames = [pd.read_parquet(f) for f in INPUTS]
pools = sorted({split_id(i)[1] for df in frames for i in df["index"]})

if os.path.exists(POOL_MAP):
    pool_map = pd.read_csv(POOL_MAP, sep="\t", dtype=str).set_index("pool")["code"].to_dict()
    missing = [p for p in pools if p not in pool_map]
    if missing:
        raise ValueError(f"pools not in {POOL_MAP}: {missing[:5]}")
else:
    codes = [f"P{i:03d}" for i in range(1, len(pools) + 1)]
    order = np.random.default_rng().permutation(len(pools))
    pool_map = {pool: codes[j] for pool, j in zip(pools, order)}
    pd.DataFrame({"pool": list(pool_map), "code": list(pool_map.values())}).to_csv(
        POOL_MAP, sep="\t", index=False
    )

for df, out in zip(frames, OUTPUTS):
    barcode_pool = df["index"].map(split_id)
    n_ids = df["index"].nunique()
    df["index"] = [f"{b}_{pool_map[p]}" for b, p in barcode_pool]
    assert df["index"].nunique() == n_ids
    df.to_parquet(out, index=False)
    print(f"wrote {out}: {len(df)} cells")
