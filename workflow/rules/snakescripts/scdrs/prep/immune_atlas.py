# script to prepare PBMCpedia data  for scDRS
import scanpy as sc
import scdrs
import pandas as pd
import yaml
import numpy as np
from pathlib import Path

INPUT = snakemake.input
OUTPUT = snakemake.output

H5AD = INPUT['h5ad']
# GS = INPUT['gs']
SAMPLE_COVAR = INPUT['sample_covar']
CONFIG = INPUT['config']

OUT_COV = OUTPUT['cov']
OUT_H5AD = OUTPUT['prep_h5ad']

# Interactive test
# H5AD = "resources/scdrs/h5ad/immune_atlas.h5ad"   
# SAMPLE_COVAR = "resources/scdrs/sample_covar/immune_atlas.sample_covar.csv"
# OUT_COV = "resources/scdrs/cov/immune_atlas.cov.tsv"
# CONFIG = "resources/scdrs/config/immune_atlas.yaml"
# OUT_H5AD = "resources/scdrs/h5ad/immune_atlas.prep.h5ad"

# Load config
with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)

# Load data
adata = sc.read(H5AD)

# count number of obs (row) per cell type, output as dictionary
d_count = adata.obs.groupby('AIFI_L2').size().to_dict()

# sample targeting 10,000 cells per cell type
d_count_sample = {k: min(v, config['max_cells_per_cell_type']) for k, v in d_count.items()}

# sample index of each cell type up to target N
sampled_index = adata.obs.groupby('AIFI_L2', group_keys=False).apply(
        lambda x: x.sample(n=d_count_sample[x.name], random_state=1234)
    ).index.tolist()

# prepare covariate
covariates = ['sample.subjectAgeAtDraw', 'sex', 'n_genes', 'disease', 'cohort.cohortGuid']
df_cov_sample = adata.obs.loc[sampled_index, covariates].copy()

# add constant = 1 column in front
df_cov_sample.insert(0, 'cons', 1)

# move index to column and change name to index
df_cov_sample = df_cov_sample.reset_index()
df_cov_sample = df_cov_sample.rename(columns={'barcodes': 'index'})

# sample adata object
adata_sample = adata[sampled_index, :].raw.to_adata()

# check if count per cell type is as expected
assert adata_sample.obs['AIFI_L2'].value_counts().to_dict() == d_count_sample

# check if sampled raw data matrix contains raw counts
assert np.all(np.mod(adata_sample.X.data, 1) == 0)

# write sampled covariate file
Path(OUT_COV).parent.mkdir(parents=True, exist_ok=True)
df_cov_sample.to_csv(OUT_COV, sep="\t", index=False)

# rwrite sampled h5ad file
Path(OUT_H5AD).parent.mkdir(parents=True, exist_ok=True)
adata_sample.write_h5ad(OUT_H5AD)