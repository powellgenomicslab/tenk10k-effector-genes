import scanpy as sc
import scdrs
import pandas as pd
import yaml
from pathlib import Path
import pickle
import numpy as np

INPUT = snakemake.input
OUTPUT = snakemake.output

H5AD = INPUT['prep_h5ad']
CONFIG = INPUT['config']
COV = INPUT['cov'] 
OUT = OUTPUT[0]

# Interactive test
# STUDY = "immune_atlas"
# H5AD = f"resources/scdrs/h5ad/{STUDY}.prep.h5ad"   
# COV = f"resources/scdrs/cov/{STUDY}.cov.tsv"
# CONFIG = f"resources/scdrs/config/{STUDY}.yaml"
# OUT = f"resources/scdrs/cov/{STUDY}.reg.h5ad.pkl"

with open(CONFIG, 'r') as f:
    config = yaml.safe_load(f)

df_cov = pd.read_csv(COV, sep="\t", index_col=0)

df_cov = pd.get_dummies(df_cov, drop_first=True).astype(float)


adata = scdrs.util.load_h5ad(
    h5ad_file=H5AD,
    flag_filter_data=config['filter_data'],
    flag_raw_count=config['raw_count']
)


# preprocess data to
# (1) regress out from covariates
# (2) group genes into bins by mean and variance
scdrs.pp.preprocess(adata, cov=df_cov,
                    n_mean_bin=config['pp_n_mean_bin'],
                    n_var_bin=config['pp_n_var_bin'],
                    copy=False)

# rewrite adata object post preprocessing
with open(OUT, 'wb') as f:
    pickle.dump(adata, f)

# Path(OUT_H5AD).parent.mkdir(parents=True, exist_ok=True)
# adata.write_h5ad(OUT_H5AD)


