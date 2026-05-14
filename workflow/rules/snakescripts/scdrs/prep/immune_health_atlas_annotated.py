# script to prepare Immune Health Atlas for scDRS
# (Annotated with TenK10K cell types using Celltypist by Peter Allen)

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
# H5AD = "resources/scdrs/h5ad/immune_health_atlas_annotated.h5ad"   
# SAMPLE_COVAR = "resources/scdrs/sample_covar/immune_health_atlas_annotated.sample_covar.csv"
# OUT_COV = "resources/scdrs/cov/immune_health_atlas_annotated.cov.tsv"
# CONFIG = "resources/scdrs/config/immune_health_atlas_annotated.yaml"
# OUT_H5AD = "resources/scdrs/h5ad/immune_health_atlas_annotated.prep.h5ad"

# Load config
with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)

# Load data
adata = sc.read(H5AD)

# Check distribution of conf_score across predicted labels
# pd.set_option('display.max_columns', None)
# pd.set_option('display.width', None)
# print("Distribution of conf_score across predicted labels:")
# print(adata.obs.groupby('predicted_labels')['conf_score'].describe())

# filter based on conf score >= 0.7
adata_filtered = adata[adata.obs['conf_score'] >= 0.8]
print(adata_filtered.obs.groupby('predicted_labels')['conf_score'].describe())

d_count = adata_filtered.obs.groupby('predicted_labels').size().to_dict()

# sample targeting 10,000 cells per cell type
d_count_sample = {k: min(v, config['max_cells_per_cell_type']) for k, v in d_count.items()}

# sample index of each cell type up to target N
sampled_index = adata_filtered.obs.groupby('predicted_labels', group_keys=False).apply(
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
# use the raw counts stored in layers.counts

adata_sample = sc.AnnData(adata_filtered[sampled_index, :].layers['counts'])

# copy all attributes in adata_filtered to adata_sample
for attr in ['uns', 'obs', 'obsm', 'obsp', 'var']:
    setattr(adata_sample, attr, getattr(adata_filtered[sampled_index, :], attr).copy())
# check data
assert adata_sample.obs['predicted_labels'].value_counts().to_dict() == d_count_sample

# check if sampled raw data matrix contains raw counts
assert np.all(np.mod(adata_sample.X.data, 1) == 0)

# write sampled covariate file
Path(OUT_COV).parent.mkdir(parents=True, exist_ok=True)
df_cov_sample.to_csv(OUT_COV, sep="\t", index=False)

# rewrite sampled h5ad file
Path(OUT_H5AD).parent.mkdir(parents=True, exist_ok=True)
adata_sample.write_h5ad(OUT_H5AD)