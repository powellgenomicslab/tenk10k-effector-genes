"""Per-cell scDeepID scores of the B cells sampled for scDRS, for the data bundle.

Follows scripts/4-polygenic/scDeepID/scDRS_scDeepID_BCELL.py (main branch). The
pathway scores of all B cells are normalised and scaled together, then restricted
to the B cells sampled for scDRS. Written per cell: the scDeepID UMAP coordinates
and cell type, and the scores of the pathways shown in Extended Data Fig 7, with the
anonymised cell IDs of the bundled scDRS files. Nothing else is shared.

Before writing, the one-sided t-test (SLE-enriched cells, scDRS Monte Carlo
P < 0.05, against the rest) is repeated for those pathways and checked against
the pathway test behind Supp Table 10, so that the cells come from the same
scDeepID run.

Run by rule scdeepid_bcell_extract in workflow/rules/scdrs.smk.
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

B_CELL_TYPES = ["B_naive", "B_intermediate", "B_memory"]
TRAIT = "sle"
PATHWAYS = list(snakemake.params.pathways)

# Pool codes must match the bundled scDRS files, so the existing map is required.
POOL_MAP = snakemake.params.pool_map
if not os.path.exists(POOL_MAP):
    raise FileNotFoundError(f"{POOL_MAP} is needed to match the bundled scDRS cell IDs")
pool_map = pd.read_csv(POOL_MAP, sep="\t", dtype=str).set_index("pool")["code"].to_dict()

score = sc.read_h5ad(snakemake.input.pathway)
sc.pp.normalize_total(score)
sc.pp.scale(score, max_value=10)
umap = sc.read_h5ad(snakemake.input.umap, backed="r")
assert (umap.obs_names == score.obs_names).all()
score.obsm["X_umap"] = np.asarray(umap.obsm["X_umap"])

drs = pd.read_parquet(snakemake.input.cell_score, columns=["index", "cell_type"])
mcp = pd.read_parquet(snakemake.input.cell_mcp, columns=["index", TRAIT]).set_index("index")
cells = drs.loc[drs["cell_type"].isin(B_CELL_TYPES), "index"]
missing = pd.Index(cells).difference(score.obs_names)
if len(missing):
    raise ValueError(f"{len(missing)} sampled B cells lack scDeepID scores, e.g. {missing[:3].tolist()}")
b = score[cells.to_numpy()].copy()
X = np.asarray(b[:, PATHWAYS].X)
enriched = mcp.loc[b.obs_names, TRAIT].to_numpy() < 0.05

reference = pd.read_csv(snakemake.input.reference).set_index("pathway").loc[PATHWAYS]
p = stats.ttest_ind(X[enriched], X[~enriched], alternative="greater").pvalue
mismatch = np.abs(np.log10(p) / np.log10(reference["p_value"].to_numpy()) - 1)
if mismatch.max() > 1e-6:
    raise ValueError(f"P values differ from {snakemake.input.reference}: is this the same scDeepID run?\n"
                     f"{pd.DataFrame({'here': p, 'reference': reference['p_value']})}")


def anonymise(cell_id):
    barcode, _, pool = cell_id.partition("_")
    return f"{barcode}_{pool_map[pool]}"


out = pd.DataFrame({
    "index": [anonymise(i) for i in b.obs_names],
    "umap_1": b.obsm["X_umap"][:, 0],
    "umap_2": b.obsm["X_umap"][:, 1],
    "scdeepid_cell_type": b.obs["celltype"].astype(str).to_numpy(),
})
out = pd.concat([out, pd.DataFrame(X, columns=PATHWAYS)], axis=1)
assert out["index"].is_unique
out.to_parquet(snakemake.output.cells, index=False)
print(f"{b.n_obs} cells, {enriched.sum()} SLE-enriched; P values match the reference")
