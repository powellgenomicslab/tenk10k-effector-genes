import scdrs
import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np
import yaml
import pickle

ADATA_PKL = snakemake.input.reg_pkl
CONFIG = snakemake.input.config
CELL_SCORE = snakemake.input.cell_score
OUTPUT = snakemake.output[0]
PARAMS = snakemake.params
INPUT_ANNOT  = snakemake.input.annot
ANNOT = snakemake.wildcards.annot

# Load config
# Interactive test
# ADATA_PKL = "resources/scdrs/h5ad/immune_atlas.reg.h5ad.pkl" 
# CELL_SCORE = "results/scdrs/cell_score/immune_atlas/crohns.cell_score.tsv.parquet.gz"
# INPUT_ANNOT = "results/scdrs/annot_sample/immune_atlas.scpred.format.csv"
# ANNOT = "scpred"
# CONFIG = "resources/scdrs/config/immune_atlas.yaml"

with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)

# Load pickled adata to get the cell type annotation
with open(ADATA_PKL, 'rb') as f:
    adata = pickle.load(f)

# compute neighborhood graph if not exist
if 'connectivities' not in adata.obsp.keys():
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)


df_celltype = pd.read_csv(INPUT_ANNOT).set_index('cell_id')

# assign cell type annotation to adata
cell_type = df_celltype.loc[adata.obs_names, 'cell_type']
adata.obs['custom_annot'] = df_celltype.loc[adata.obs_names, 'cell_type']

df_score = pd.read_parquet(CELL_SCORE)

# compute cell type stats
df_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["custom_annot"],
)["custom_annot"] 

Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)
df_stats \
    .reset_index() \
    .to_csv(OUTPUT, sep='\t', index=False)