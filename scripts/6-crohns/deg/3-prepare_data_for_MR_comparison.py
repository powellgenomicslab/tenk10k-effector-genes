# prepare adata colon for MR comparison 

# qsub -I -q normal -P fy54 -l ncpus=4,storage=gdata/fy54+gdata/ei56,mem=50GB -l jobfs=100GB
# conda activate scanpy

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt

import os

sc.set_figure_params(dpi=100, dpi_save=300, color_map="viridis_r")
sc.settings.verbosity = 0
sc.logging.print_header()

# This doesn't work 
# print(plt.rcParams["font.sans-serif"])
# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams["font.sans-serif"] = ['Helvetica']

location = "colon"

# Set directories 
home_dir= "/g/data/ei56/rt3501"
data_dir = f"{home_dir}/data/crohnsdata/{location}"
project_dir = f"{home_dir}/crohns_vignette"
figs = f"{project_dir}/preparing_data/figs"

adata = sc.read_h5ad(f"{data_dir}/cd_colon_immune.h5ad")

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata, 
        color="Celltype", 
        title = "Original Annotation",
        frameon=False, 
        vmax="p99", 
        legend_fontsize=10,
        legend_fontoutline=2,
        s=20)
plt.savefig(f"{figs}/crohns_umap_original_cell_annotation.png", bbox_inches='tight', dpi=300)


# change 'normal' to healthy
# Create a new column called condition status which is the same as Type, but with more words
# Update the disease column to reflect this, as the current disease column is mislabelled
adata.obs['condition_status'] = ''
adata.obs.loc[adata.obs['Type'] == 'Heal', 'condition_status'] = 'Healthy Colon' 
adata.obs.loc[adata.obs['Type'] == 'Infl', 'condition_status'] = 'Inflamed Crohns Disease Colon' 
adata.obs.loc[adata.obs['Type'] == 'NonI', 'condition_status'] = 'Non-inflamed Crohns Disease Colon' 

types = ['Healthy Colon', 'Non-inflamed Crohns Disease Colon', 'Inflamed Crohns Disease Colon']

adata.obs['disease']

# Categories (2, object): ['Crohn disease', 'normal']
#fix the disease col
adata.obs['disease_corrected'] = ''
adata.obs.loc[adata.obs['Type'] == 'Heal', 'disease_corrected'] = 'Healthy' 
adata.obs.loc[adata.obs['Type'] == 'Infl', 'disease_corrected'] = 'Crohns Disease' 
adata.obs.loc[adata.obs['Type'] == 'NonI', 'disease_corrected'] = 'Crohns Disease' 
adata.obs['disease_corrected'] = adata.obs['disease_corrected'].astype('category')

# note the deg results cell types are not used, i took the adata colon cell types 
major_cell_types = pd.read_csv(f"{data_dir}/adata_colon_celltype_groups.csv")

major_cell_types.columns = ["original", "major_cell_type"]

# used AI for this dictionary 
id_to_symbol = dict(zip(major_cell_types["original"], major_cell_types["major_cell_type"]))

# this gets the original cell types and maps them to the major cell types
mapped_cts = [id_to_symbol.get(original, original) for original in adata.obs['Celltype']]

# create a new column in the adata object for major cell types
adata.obs['major_cell_type'] = mapped_cts
adata.obs['major_cell_type'] = adata.obs['major_cell_type'].astype('category')

###############################################################################################
# add genes symbols 
# get gene names from features file of raw data. you can also just string split the features col in adata (formatted as ENSID_genesymbol) to pull out the gene symbol but idk python
features = pd.read_csv(f"{data_dir}/raw_data/imm/features_original.tsv", header=None, sep = "\t")
features.columns = ["ensembl_id", "gene_symbol"]

# used AI for this dictionary 
# Create a mapping dictionary from Ensembl ID to gene symbol
id_to_symbol = dict(zip(features["ensembl_id"], features["gene_symbol"]))

# Map the Ensembl IDs in adata.var_names to gene symbols
mapped_symbols = [id_to_symbol.get(ensembl_id, ensembl_id) for ensembl_id in adata.var_names]

# Update adata.var_names with gene symbols
adata.var['gene_symbols']= mapped_symbols

################################################################################################################################
# save the object with major cell types and gene_symbol columns 
filename = f"{home_dir}/cd_colon_immune_major_cell_types.h5ad"
adata.write_h5ad(filename)
################################################################################################################################
