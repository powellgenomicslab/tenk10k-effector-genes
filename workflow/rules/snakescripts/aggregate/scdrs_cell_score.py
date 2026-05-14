import scdrs
import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np
import yaml
import pickle

ADATA = snakemake.input.prep_h5ad
CONFIG = snakemake.input.config
CELL_SCORE_FILES = snakemake.input.cell_score
OUTPUT = snakemake.output

# interactive
# STUDY = "immune_atlas"
# ADATA_PKL = f"resources/scdrs/h5ad/{STUDY}.reg.h5ad.pkl"
# ADATA = f"resources/scdrs/h5ad/{STUDY}.prep.h5ad"
# CONFIG = f"resources/scdrs/config/{STUDY}.yaml"
# CELL_SCORE_FILES = list(Path(f"results/scdrs/cell_score/{STUDY}").glob("*.cell_score.tsv.parquet.gz"))
# OUTPUT = f"results/aggregate/{STUDY}/scdrs.cell_score.tsv.parquet.gz"

# combine all scores into a single DataFrame
def combine_scores(col):
    """Combine scores from multiple cell score files into a single DataFrame."""
    df = pd.concat([
        pd.read_parquet(f, columns=[col]).rename(
            columns={col: Path(f).name.removesuffix('.cell_score.tsv.parquet.gz')}
        )
        for f in CELL_SCORE_FILES
    ], axis=1)
    return df

df_scores = combine_scores('norm_score')

df_mcp = combine_scores('mc_pval')

with open(CONFIG, "r") as f:
    config = yaml.safe_load(f)

adata = sc.read_h5ad(ADATA)

# preprocessing
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, config['regress_out_cols'])
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver="arpack", n_comps=30)

# perform Harmony integration to remove technical batch effects between pools
sc.external.pp.harmony_integrate(adata, config['harmony_batch_col'])

# Re-run UMAP on the Harmony principal components
adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)

# Get annotation (cell type) and recomputed UMAP coordinates

df_annot = adata \
    .obs[[config['group_stats_col']]] \
    .rename(columns={config['group_stats_col']: 'cell_type'})
df_umap = pd.DataFrame(adata.obsm['X_umap'],
                    index=adata.obs_names,
                    columns=['umap_1', 'umap_2'])

# Combine scores with annotations & umap & save
df_merged = pd.concat(
    [df_annot, df_umap, df_scores],
    axis=1,
    join='inner'
).reset_index()

df_merged.to_parquet(OUTPUT['scores'], index=False, compression='gzip')

# Combine cell-type mc P-value with annotations & umap & save
df_merged = pd.concat(
    [df_annot, df_umap, df_mcp],
    axis=1,
    join='inner'
).reset_index()
df_merged.to_parquet(OUTPUT['mcp'], index=False, compression='gzip')