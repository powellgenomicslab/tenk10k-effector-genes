# Purpose: Prepare supplementary fig - UMAP of original annotation and major cell types shared with TenK10K 

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
from matplotlib import font_manager
import os

sc.set_figure_params(dpi=100, dpi_save=300, color_map="viridis")
sc.settings.verbosity = 0
sc.logging.print_header()

# This doesn't work 
# print(plt.rcParams["font.sans-serif"])
# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams["font.sans-serif"] = ['Helvetica']
# plt.rcParams["font.family"] = "Helvetica"

# Set directories 
data_dir = "resources/crohns_case_study/deg"
figs = "resources/crohns_case_study/figures/deg"

# get Helvetica 
font_path = '/g/data/ei56/rt3501/miniforge3/envs/scanpy/fonts/Helvetica.ttf'  # Your font path goes here
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()

adata = sc.read_h5ad(f"{data_dir}/cd_colon_immune_major_cell_types.h5ad")

################################################################################################################################################################################################################################
# Plot Cell Annotations 
################################################################################################################################################################################################################################

np.random.seed(525)

alpha_var = 1
s_var = 3

# colours = ['Celltype', 'major_cell_type']
# titles = ["Original Annotation","Major Cell Types"]

colours = ['Celltype']
titles = ["Original Annotation"]

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata, color=colours,
        frameon=False, 
        ncols = 2,
        title = titles,
        alpha = alpha_var,
        vmax="p99", 
        legend_fontsize=5,
        legend_loc='right margin',
        #legend_loc='lower right',
        # add_outline=True,        # Adds a nice 'finished' look to clusters
        # outline_width=(0.5, 0.05),
        legend_fontoutline=2,
        s=s_var)

# plt.savefig(f"{figs}/original_and_major_annotation_all_celltypes_alpha_0.5.png", bbox_inches='tight', dpi=300)
plt.savefig(f"{figs}/original_annotation_all_celltypes_alpha{alpha_var}_s{s_var}_add_outline.png", bbox_inches='tight', dpi=300)
plt.savefig(f"{figs}/original_annotation_all_celltypes_alpha{alpha_var}_s{s_var}.png", bbox_inches='tight', dpi=300)

# Grey out the major cell types that are not shared in the TenK10K dataset 
groups = ['CD4 T', 'CD8 T', 'NK', 'Unconventional T', 'B', 'Plasma B', 'Monocyte', 'Dendritic']

# Because of 'groups' nothing shows up under original annotation
colours = ['major_cell_type']
titles = ["Major Cell Types"]

with rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata, 
        color=colours, 
        groups=groups, # filter cell types 
        frameon=False, 
        title = titles,
        ncols = 2,
        vmax="p99", 
        alpha = alpha_var,
        legend_fontsize=5,
        legend_loc='right margin',
        # add_outline=True,        # Adds a nice 'finished' look to clusters
        # outline_width=(0.5, 0.05),
        legend_fontoutline=2,
        s=s_var)

# plt.savefig(f"{figs}/original_and_major_annotation_filtered_celltypes_alpha{alpha_var}_s{s_var}_add_outline.png", bbox_inches='tight', dpi=300)

plt.savefig(f"{figs}/original_and_major_annotation_filtered_celltypes_alpha{alpha_var}_s{s_var}.png", bbox_inches='tight', dpi=300)
