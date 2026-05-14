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

# Load config
with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)

CELL_TYPE_COL = config['group_stats_col']
# Load pickled adata to get the cell type annotation
with open(ADATA_PKL, 'rb') as f:
    cell_type = pickle.load(f).obs[CELL_TYPE_COL]

cols = ['raw_score', 'norm_score', 'mc_pval', 'pval', 'nlog10_pval', 'zscore']

df = pd.merge(
    pd.read_parquet(CELL_SCORE).loc[:, cols],
    cell_type,
    left_index=True,
    right_index=True
).reset_index()

# count N cells in top 95th percentile and 99th percentile per cell type
percentiles = PARAMS['percentiles']

counts = {}
for percentile in percentiles:
    percentile = int(percentile)
    threshold = np.percentile(df['norm_score'], percentile)
    filtered_df = df[df['norm_score'] >= threshold]
    counts[percentile] = filtered_df \
        .groupby(CELL_TYPE_COL) \
        .size() \
        .rename(f'top_{percentile}th_pct_count') \
        .rename_axis('cell_type')

# Combine counts into a single dataframe
counts_df = pd.concat(counts.values(), axis=1).reset_index()

# Save counts to a file
Path(str(OUTPUT)).parent.mkdir(parents=True, exist_ok=True)

counts_df.to_csv(str(OUTPUT), sep="\t", index=False)    
