import scanpy as sc
import scdrs
import pandas as pd
import yaml
from pathlib import Path

INPUT = snakemake.input
OUTPUT = snakemake.output

H5AD = INPUT['h5ad']
SAMPLE_COVAR = INPUT['sample_covar']

OUT_COV = OUTPUT['cov']

# Interactive test
# H5AD = "resources/scdrs/h5ad/tenk10k_phase1.h5ad"   
# SAMPLE_COVAR = "resources/scdrs/sample_covar/tenk10k_phase1.sample_covar.csv"
# OUT_COV = "resources/scdrs/cov/tenk10k_phase1.cov.tsv"
# OUT_H5AD = "resources/scdrs/h5ad/tenk10k_phase1.regressed.h5ad"

col_sample_covar = ['sample_id', 'sex', 'age']
df_sample = pd.read_csv(SAMPLE_COVAR).loc[:, col_sample_covar]

# Load data
adata = sc.read_h5ad(H5AD)

# get the barcode level metadata 
metadata = adata.obs.loc[:,["wg2_scpred_prediction", "cpg_id", "n_genes_by_counts", "cohort"]]

df_cov = pd.merge(
    metadata.reset_index(),
    df_sample, 
    left_on='cpg_id', 
    right_on='sample_id'
)

df_cov['cons'] = 1
df_cov['sex'] = df_cov['sex'].astype('int') - 1
df_cov['age'] = df_cov['age'].astype('int')
df_cov['cohort'] = [1 if x == 'TOB' else 0 for x in df_cov['cohort']]

select_cols = ['index', 'cons', 'sex', 'age', 'n_genes_by_counts', 'cohort']

df_cov = df_cov.loc[:,select_cols]

# write
Path(OUT_COV).parent.mkdir(parents=True, exist_ok=True)
df_cov.to_csv(OUT_COV, sep="\t", index=False)


# df_cov = df_cov.set_index('index')
# # load geneset, convert homologs and overlap gene names to adata.var_names
# dict_gs = scdrs.util.load_gs(
#     GS,
#     src_species=config['gs_src_species'],
#     dst_species=config['gs_dst_species'],
#     to_intersect=adata.var_names,
# )

# # preprocess data to
# # (1) regress out from covariates
# # (2) group genes into bins by mean and variance
# scdrs.pp.preprocess(adata, cov=df_cov,
#                     n_mean_bin=config['pp_n_mean_bin'],
#                     n_var_bin=config['pp_n_var_bin'],
#                     copy=False)

# # rewrite adata object post preprocessing
# Path(OUT_H5AD).parent.mkdir(parents=True, exist_ok=True)
# adata.write_h5ad(OUT_H5AD)


