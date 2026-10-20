"""Compare scDeepID cell function scores between SLE-enriched and other B cells.

Ported from sections 2-7 of the scDeepID revision notebook,
revision_cleanup_v3_github.ipynb (copied, without outputs, next to this script);
its figures are redrawn by scripts/6-polygenic/scDeepID/. For each cell function
(GO biological process), the scores of all B cells are normalised and scaled
together, restricted to the B cells sampled for scDRS, and compared between
SLE-enriched cells (scDRS Monte Carlo P < 0.05) and the rest with a one-sided
t-test, with Benjamini-Hochberg FDR. Each cell function is then flagged with the
GO categories (children of "biological process") that contain it, keeping the
categories with at least MIN_PATHWAYS cell functions at FDR < 0.05.

Writes pathway_test_full.csv (Supp Table 10, Supp Fig 15, Extended Data Fig 7).
Run by rule scdeepid_pathway_test in workflow/rules/scdrs.smk.
"""

import re

import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

B_CELL_TYPES = ["B_naive", "B_intermediate", "B_memory"]
BP_ROOT = "GO:0008150"
MIN_PATHWAYS = 5

# Scores and scDRS ---------------------------------------------------------------
pathway_attr = sc.read(snakemake.input.pathway)
sc.pp.normalize_total(pathway_attr)
sc.pp.scale(pathway_attr, max_value=10)

keep_cols = ["index", "cell_type", "sle"]
df_score = pd.read_parquet(snakemake.input.cell_score, columns=keep_cols)
df_mcp = pd.read_parquet(snakemake.input.cell_mcp, columns=keep_cols)
df_score = df_score[df_score["cell_type"].isin(B_CELL_TYPES)]
df_mcp = df_mcp[df_mcp["cell_type"].isin(B_CELL_TYPES)].set_index("index")

common_cells = sorted(set(pathway_attr.obs_names) & set(df_score["index"]))
subset = pathway_attr[common_cells].copy()
sle_status = np.where(df_mcp.loc[subset.obs_names, "sle"].to_numpy() < 0.05, "SLE", "Non-SLE")

# Per-pathway t-test and FDR -----------------------------------------------------
p_values, sle_means, non_sle_means = [], [], []
for pathway in subset.var_names:
    sle_data = subset[sle_status == "SLE", pathway].X.flatten()
    non_sle_data = subset[sle_status == "Non-SLE", pathway].X.flatten()
    _, p_val = stats.ttest_ind(sle_data, non_sle_data, alternative="greater")
    p_values.append(p_val)
    sle_means.append(sle_data.mean())
    non_sle_means.append(non_sle_data.mean())

_, fdr_values, _, _ = multipletests(p_values, method="fdr_bh")
results = pd.DataFrame({
    "pathway": list(subset.var_names),
    "neg_log_p": -np.log10(np.array(p_values)),
    "p_value": p_values,
    "fdr": fdr_values,
    "sle_mean": sle_means,
    "non_sle_mean": non_sle_means,
}).sort_values("neg_log_p", ascending=False).reset_index(drop=True)


# GO categories --------------------------------------------------------------------
def parse_obo(path):
    terms, cur = {}, None
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line == "[Term]":
                cur = {"parents": [], "namespace": None, "name": None, "obsolete": False}
                continue
            if line.startswith("[") and line != "[Term]":
                if cur and cur.get("id"):
                    terms[cur["id"]] = cur
                cur = None
                continue
            if cur is None:
                continue
            if line == "":
                if cur.get("id"):
                    terms[cur["id"]] = cur
                cur = None
                continue
            if line.startswith("id: "):
                cur["id"] = line[4:]
            elif line.startswith("name: "):
                cur["name"] = line[6:]
            elif line.startswith("namespace: "):
                cur["namespace"] = line[11:]
            elif line.startswith("is_a: "):
                cur["parents"].append(line[6:].split(" ! ")[0])
            elif line.startswith("relationship: part_of "):
                cur["parents"].append(line[len("relationship: part_of "):].split(" ! ")[0])
            elif line == "is_obsolete: true":
                cur["obsolete"] = True
    if cur and cur.get("id"):
        terms[cur["id"]] = cur
    return terms


terms = parse_obo(snakemake.input.go_obo)
G = nx.DiGraph()
for tid, t in terms.items():
    if t.get("obsolete"):
        continue
    for p in t["parents"]:
        G.add_edge(p, tid)


def msigdb_name(name):
    return "GOBP_" + re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").upper()


def slugify(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def bp_descendants(root):
    d = nx.descendants(G, root) | {root}
    return {x for x in d
            if terms.get(x, {}).get("namespace") == "biological_process"
            and not terms.get(x, {}).get("obsolete")
            and terms.get(x, {}).get("name")}


sig = set(results.loc[results["fdr"] < 0.05, "pathway"])
categories = []
for cid in sorted(G.successors(BP_ROOT)):
    info = terms.get(cid, {})
    if info.get("namespace") != "biological_process" or info.get("obsolete") or not info.get("name"):
        continue
    members = {msigdb_name(terms[d]["name"]) for d in bp_descendants(cid)}
    if len(sig & members) >= MIN_PATHWAYS:
        categories.append((len(sig & members), info["name"], members))
categories.sort(key=lambda c: -c[0])

for n_sig, name, members in categories:
    results[f"cat_{slugify(name)}"] = results["pathway"].isin(members)
    print(f"{name}: {n_sig} cell functions at FDR < 0.05")

results.to_csv(snakemake.output.test, index=False)
print(f"{len(common_cells)} cells, {int((sle_status == 'SLE').sum())} SLE-enriched; "
      f"{len(results)} cell functions, {len(sig)} at FDR < 0.05")
