# umap of the original cell annotation: all and filtered to MR matched 
# umap of the major cell type annotation: all and filtered to MR matched

# get the umap of the genes of interest, by disease category only as well as two subtypes: all and filtered to MR matched

# use the font helvetica 

# qsub -I -q normal -P fy54 -l ncpus=4,storage=gdata/fy54+gdata/ei56,mem=50GB -l jobfs=100GB

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
from matplotlib import font_manager
import os

sc.set_figure_params(dpi=100, dpi_save=300, color_map="viridis_r")
sc.settings.verbosity = 0
sc.logging.print_header()

# This doesn't work 
# print(plt.rcParams["font.sans-serif"])
# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams["font.sans-serif"] = ['Helvetica']
# plt.rcParams["font.family"] = "Helvetica"

location = "colon"

# Set directories 
home_dir= "/g/data/fy54/analysis/tenk10k-causal"
data_dir = f"{home_dir}/resources/crohns_case_study"
figs = "/g/data/ei56/rt3501/crohns_vignette/deg/figs_main_or_supplementary"

# get Helvetica 
font_path = '/g/data/ei56/rt3501/miniforge3/envs/scanpy/fonts/Helvetica.ttf'  # Your font path goes here
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()

adata = sc.read_h5ad(f"{data_dir}/deg/cd_colon_immune_major_cell_types.h5ad")

################################################################################################################################################################################################################################
# Plot Cell Annotations 
################################################################################################################################################################################################################################

np.random.seed(525)
# random_indices = np.random.permutation(list(range(adata.shape[0])))

# Plot ALL cells with original and major cell type annotations  
alpha_var = 1
s_var = 3

colours = ['major_cell_type', 'Celltype']
titles = ["Major Cell Types", "Original Annotation"]

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata[random_indices, :], color=colours,
        frameon=False, 
        ncols = 2,
        title = titles,
        alpha = alpha_var,
        vmax="p99", 
        legend_fontsize=5,
        legend_loc='right margin',
        #legend_loc='lower right',
        legend_fontoutline=2,
        s=s_var)
#plt.savefig(f"{figs}/original_and_major_annotation_all_celltypes_alpha_0.5.png", bbox_inches='tight', dpi=300)
plt.savefig(f"/g/data/ei56/rt3501/original_and_major_annotation_all_celltypes_alpha{alpha_var}_s{s_var}.png", bbox_inches='tight', dpi=300)

# You can't make a UMAP on filtered data so ignore this 
# adata_filtered = adata[~adata.obs['major_cell_type'].isin(['Immune Cycling cells', 'Macrophage', 'Mast', 'Intraepithelial NK-like']), :]

# Plot only cells that match MR cell types (using major cell type annotations) with original and major cell type annotations. 

colours = ['major_cell_type', 'Celltype']
titles = ["Major Cell Types", "Original Annotation"]

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata[~adata.obs['major_cell_type'].isin(['Immune Cycling cells', 'Macrophage', 'Mast', 'Intraepithelial NK-like'])], color=colours,
        frameon=False, 
        title = titles,
        ncols = 2,
        vmax="p99", 
        alpha = alpha_var,
        legend_fontsize=5,
        #legend_loc='right margin',
        legend_loc='right margin',
        legend_fontoutline=2,
        s=s_var)
#plt.savefig(f"{figs}/original_and_major_annotation_filtered_celltypes_alpha{alpha_var}_s{s_var}.png", bbox_inches='tight', dpi=300)
# Inode error saving it to correct directory so have to export to home dir
plt.savefig(f"/g/data/ei56/rt3501/original_and_major_annotation_filtered_celltypes_alpha{alpha_var}_s{s_var}.png", bbox_inches='tight', dpi=300)

# Plot only cells that match MR cell types (using major cell type annotations) major cell type annotations only. 
alpha_var = 1
s_var = 3
with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata[~adata.obs['major_cell_type'].isin(['Immune Cycling cells', 'Macrophage', 'Mast', 'Intraepithelial NK-like'])], color=['Celltype'],
        frameon=False, 
        title = ["Original"],
        alpha = alpha_var,
        ncols = 1,
        vmax="p99", 
        legend_fontsize=5,
        #legend_loc='right',
        legend_fontoutline=2,
        s=s_var)
plt.savefig("/g/data/ei56/rt3501/original_annotation_filtered_cell_types_s3.png", bbox_inches='tight', dpi=300)

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata, color=['Celltype'],
        frameon=False, 
        title = ["Original"],
        alpha = alpha_var,
        ncols = 1,
        vmax="p99", 
        legend_fontsize=5,
        #legend_loc='right',
        legend_fontoutline=2,
        s=s_var)
plt.savefig("/g/data/ei56/rt3501/original_annotation_all_cell_types_s3.png", bbox_inches='tight', dpi=300)

# Inode error saving it to correct directory below
# plt.savefig(f"{figs}/major_cell_type_annotation_original_and_major_filtered_cell_types.png", bbox_inches='tight', dpi=300)
