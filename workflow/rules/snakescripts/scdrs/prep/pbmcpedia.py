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
# H5AD = "resources/scdrs/h5ad/pbmcepedia.h5ad"   
# SAMPLE_COVAR = "resources/scdrs/sample_covar/pbmcpedia.sample_covar.csv"
# OUT_COV = "resources/scdrs/cov/pbmcpedia.cov.tsv"
# CONFIG = "resources/scdrs/config/pbmcpedia.yaml"
# OUT_H5AD = "resources/scdrs/h5ad/pbmcpedia.prep.h5ad"


# Load data
adata = sc.read(H5AD)

# chec[k dimension of adata obs
print(f"adata obs dim: {adata.obs.shape}")

# prepare covariate
covariates = ['donor_id', 'sex', 'age', 'disease', 'project']
donor_meta = adata.obs[covariates].drop_duplicates().reset_index(drop=True)

# remove missing sex / age data
missing_str = ['unknown', 'nan', 'na']
donor_meta['age'] = pd.to_numeric(donor_meta['age'], errors='coerce')
donor_meta = donor_meta[~donor_meta['sex'].isin(missing_str) & ~donor_meta['age'].isna()]
donor_meta['sex'] = [1 if x == 'female' else 0 for x in donor_meta['sex'].tolist()]


# filter adata to only donors with non-missing sex / age data


adata_sample = adata[adata.obs['donor_id'].isin(donor_meta['donor_id']), :]

donor_meta = donor_meta[donor_meta.duplicated(subset=['donor_id'])]


# discard age that contains non-digit characters or .
donor_meta = donor_meta[~donor_meta['age'].str.contains(r'\D')]

donor_meta = donor_meta[donor_meta['age'].str.match(r'^\d+$')]

# check type of each columns
dtype_dict = donor_meta.dtypes.to_dict()
donor_meta['age'] = donor_meta['age'].astype('int')

# Load config
with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)
# Set index to ensgid
adata.var.set_index(['gene_id'], inplace=True, drop=False)

# Prepare covariate
# get the barcode level metadata 
metadata = adata.obs.loc[:,["wg2_scpred_prediction", "cpg_id", "n_genes_by_counts", "cohort"]]

col_sample_covar = ['donor_id', 'sex', 'age']
df_sample = pd.read_csv(SAMPLE_COVAR).loc[:, col_sample_covar]

df_cov = pd.merge(
    metadata.reset_index(),
    df_sample, 
    left_on='cpg_id', 
    right_on='sample_id'
)

# count number of obs per cell type, output as dictionary
d_count = df_cov.groupby('wg2_scpred_prediction')['index'].count().to_dict()

# sample targeting 10,000 cells per cell type
d_count_sample = {k: min(v, config['max_cells_per_cell_type']) for k, v in d_count.items()}

# downsample covariates according to target cell per cell type
df_cov_sample = df_cov \
    .groupby('wg2_scpred_prediction', observed=True) \
    .apply(lambda x: x.sample(n=d_count_sample[x.name], random_state=1234),
              include_groups='wg2_scpred_prediction') \
    .reset_index(level=0, drop=True)

# downsample covariates & data
# df_cov = df_cov \
#     .groupby(config['sample_by_col'], observed=True) \
#     .apply(lambda x: x.sample(frac=config['sample_fraction'], random_state=1234),
#            include_groups='wg2_scpred_prediction') \
#     .reset_index(level=0, drop=True)

df_cov_sample['cons'] = 1
df_cov_sample['sex'] = df_cov_sample['sex'].astype('int') - 1
df_cov_sample['age'] = df_cov_sample['age'].astype('int')
df_cov_sample['cohort'] = [1 if x == 'TOB' else 0 for x in df_cov_sample['cohort']]

select_cols = ['index', 'cons', 'sex', 'age', 'n_genes_by_counts', 'cohort']

df_cov_sample = df_cov_sample.loc[:,select_cols]

# write
Path(OUT_COV).parent.mkdir(parents=True, exist_ok=True)
df_cov_sample.to_csv(OUT_COV, sep="\t", index=False)

# sample adata
index_sample = df_cov_sample['index'].tolist()

adata_sample = adata[index_sample, :].raw.to_adata()

# rewrite adata object
Path(OUT_H5AD).parent.mkdir(parents=True, exist_ok=True)
adata_sample.write_h5ad(OUT_H5AD)