# Purpose: Annotate adata object with major cell type annotation
# Make a corrected disease status column ('disease_corrected') due to the bug in original dataset from CELLXGENE - Type == NonI was labelled "normal" under the "disease" column 
# Make a verbose column for Type (NonI, Infl, Heal)

# qsub -I -q normal -P fy54 -l ncpus=2,storage=gdata/fy54+gdata/ei56,mem=50GB -l jobfs=100GB
# cd /g/data/fy54/rt3501/repos/tenk10k-causal/
# conda activate scanpy
# python

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt

import os

sc.set_figure_params(dpi=100, dpi_save=300, color_map="viridis")
sc.settings.verbosity = 0
sc.logging.print_header()


# Set directories 

data_dir = "resources/crohns_case_study/deg"
figs = "resources/crohns_case_study/figures/deg"

adata = sc.read_h5ad(f"{data_dir}/cd_colon_immune.h5ad")

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
# Fix the disease col
adata.obs['disease_corrected'] = ''
adata.obs.loc[adata.obs['Type'] == 'Heal', 'disease_corrected'] = 'Healthy' 
adata.obs.loc[adata.obs['Type'] == 'Infl', 'disease_corrected'] = 'Crohns Disease' 
adata.obs.loc[adata.obs['Type'] == 'NonI', 'disease_corrected'] = 'Crohns Disease' 
adata.obs['disease_corrected'] = adata.obs['disease_corrected'].astype('category')

# new major cell type annotation. print out the cell types  
celltypes = pd.DataFrame(adata.obs['Celltype'].unique(), columns=['celltypes'])
celltypes.sort_values(by = 'celltypes')

# Define major cell type conditions as a list of (condition, value) tuples
conditions = [
    (adata.obs['Celltype'].isin(['Tregs',"T cells CD4 FOSB", "T cells Naive CD4", "T cells CD4 IL17A"]), "CD4 T"),
    (adata.obs['Celltype'].isin(["T cells CD8", "T cells CD8 KLRG1"]), "CD8 T"),
    (adata.obs['Celltype'].isin(["T cells OGT"]), "Unconventional T"),
    (adata.obs['Celltype'].isin(["B cells"]), "B"),
    (adata.obs['Celltype'].isin(["Plasma cells"]), "Plasma B"),
    (adata.obs['Celltype'].isin(["ILCs", "NK cells KLRF1 CD3G-"]), "NK"),
    (adata.obs['Celltype'].isin(["Monocytes CHI3L1 CYP27A1", "Monocytes S100A8 S100A9"]), "Monocyte"),
    (adata.obs['Celltype'].isin(["DC1", "DC2 CD1D", "DC2 CD1D-", "Mature DCs"]), "Dendritic"),
    (adata.obs['Celltype'].isin(["Macrophages LYVE1", "Macrophages Metallothionein", "Macrophages CCL3 CCL4", "Macrophages"]), "Macrophage"),
    (adata.obs['Celltype'].isin(["NK-like cells ID3 ENTPD1"]), "Intra-epithelial Lymphocytes"),
    (adata.obs['Celltype'].isin(["Immune Cycling cells"]), "Cycling"),
    (adata.obs['Celltype'].isin(["B cells AICDA LRMP"]), "Germinal Centre B"),
    (adata.obs['Celltype'].isin(["Mast"]), "Mast")
]

adata.obs['major_cell_type'] = adata.obs['Celltype'].case_when(conditions)
adata.obs['major_cell_type'] = adata.obs['major_cell_type'].astype('category')

############################################################################################### 
# add genes symbols 
# get gene names from features file of raw data. you can also just string split the features col in adata (formatted as ENSID_genesymbol) to pull out the gene symbol but idk python
# features = pd.read_csv(f"{data_dir}/features_original.tsv", header=None, sep = "\t")
# features.columns = ["ensembl_id", "gene_symbol"]

# # used AI for this dictionary 
# # Create a mapping dictionary from Ensembl ID to gene symbol
# id_to_symbol = dict(zip(features["ensembl_id"], features["gene_symbol"]))

# # Map the Ensembl IDs in adata.var_names to gene symbols
# mapped_symbols = [id_to_symbol.get(ensembl_id, ensembl_id) for ensembl_id in adata.var_names]

# # Update adata.var_names with gene symbols
# adata.var['gene_symbols']= mapped_symbols

################################################################################################################################
# save the object with major cell types  
filename = f"{data_dir}/cd_colon_immune_major_cell_types.h5ad"
adata.write_h5ad(filename)
################################################################################################################################
