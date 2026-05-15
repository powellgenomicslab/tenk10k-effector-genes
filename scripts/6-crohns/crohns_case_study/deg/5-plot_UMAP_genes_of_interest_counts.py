# Purpose: Plot counts for genes of interest 

# qsub -I -q normal -P fy54 -l ncpus=2,storage=gdata/fy54+gdata/ei56,mem=32GB -l jobfs=50GB
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

# Set directories 
data_dir = "resources/crohns_case_study/deg"
figs = "resources/crohns_case_study/figures/deg"

# get Helvetica 
font_path = '/g/data/ei56/rt3501/miniforge3/envs/scanpy/fonts/Helvetica.ttf' 
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()

# Load adata (pre-annotated/not!)
adata = sc.read_h5ad(f"{data_dir}/cd_colon_immune_major_cell_types.h5ad")
# create a gene_symbols col for plotting. Feature_names are ensemblids_genesymbols.
adata.var['gene_symbols'] = adata.var['feature_name'].str.split('_').str[0]

################################################################################################################################################################################################################################
# Plot Genes of Interest
# Either by Type (adata.obs.condition_status, nCat = 3) or Disease (adata.obs.disease_corrected, nCat = 2). See prepare_data_for_MR_comparison script.
################################################################################################################################################################################################################################
# Plot with either all cells (adata) or subset to shared major cell types: 

groups = ['CD4 T', 'CD8 T', 'NK', 'Unconventional T', 'B', 'Plasma B', 'Monocyte', 'Dendritic']
# groups = ['Tregs', "T cells CD4 FOSB", "T cells Naive CD4", "T cells CD4 IL17A", 
#     "T cells CD8", "T cells CD8 KLRG1", "T cells OGT", 
#     "B cells", "Plasma cells", "ILCs", "NK cells KLRF1 CD3G-", 
#     "Monocytes CHI3L1 CYP27A1", "Monocytes S100A8 S100A9", 
#     "DC1", "DC2 CD1D", "DC2 CD1D-", "Mature DCs"]
adata_subset = adata[adata.obs['major_cell_type'].isin(groups)].copy()

diseases = ['Healthy', 'Crohns Disease']
types = ['Healthy Colon', 'Non-inflamed Crohns Disease Colon', 'Inflamed Crohns Disease Colon']

# Set the axes of the final figure = columns are the conditions, rows are the genes. 
# For diseases ncols = 2, for types, ncol =3
# adjust wspace for space between figures 

############# Decide on your genes dictionary 
# genes_of_int = {'con_pos': ["SLC39A8", "CDC42SE2", "FUBP1", "PPP1R14B", "RSBN1", "LNPEP", "IRF1", "LRRK2", "ITGA4", "IPMK", "FAM53B", "SPATS2L", "ZBTB38", "HSPD1", "LPP", "CAST", "RASGRP1", "TNFRSF9", "RBM6", "ATXN2L", "JAK2"], 
#             'con_neg': ["PSMB9", "CFL1", "RNASET2", "NCF4", "PARK7", "RPS9", "CTSW", "NDFIP1", "RHOC", "RGS14", "UBE2L3", "GPX1", "HLA-A", "PTGER4", "M6PR", "SOCS1", "LITAF", "PSTPIP1"] 
#             'dis_mr_pos': ["GPX4", "TUFM", "RSBN1", "CISD1", "TMEM258", "SMIM19", "SNN", "KDELR2", "IFITM2", "IL12RB1", "BIN1", "GPX1"], 
#             'dis_mr_neg': ["PPP3CA", "IL26", "SENP7", "ASH1L", "NFATC2IP", "SYT11", "THEMIS", "TNFRSF18", "CNN2", "RAD50", "PTPN22", "ZMIZ1", "SERBP1", "PTPRC", "PSTPIP1"], 
#  }

# d_genes = {'NCF4': diseases,
#             'ZBTB38': diseases}

# d_genes = {'AIF1': types, #MATCHED DIR
#             'PTGER4': types, #MATCHED DIR
#             'CDC42SE2': types, #MATCHED DIR
#             'GPX1': types, # MIXED MATCHED DIR
#             'TNFRSF18': types, # False DIR
#             'BIN1': types} # False DIR

# d_genes = {'ITGA4': types, #MATCHED DIR
#             'IRF1': types, #MATCHED DIR
#             'TNFRSF18': types, # UNMATCHED DIR
#             'GPX1': types, # MIXED MATCHED DIR depending on major cell type
#             'TNFRSF18': types, # MIXED MATCHED DIR depending on major cell type
#             'PSTPIP1': types} # MIXED MATCHED DIR depending on major cell type


# All T cell genes? Filter adata_subset to just T cells?

d_genes = {'NCF4': types,
            'ZBTB38': types, 
            'major_cell_type': types}

############# Set the figure columns - columns are different 'types', rows are genes 

num_col = len(types)
num_row = len(d_genes)

ncols = num_col
nrows = num_row 
fsize = 6
wspace = 0.3
fig, axs = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(ncols * fsize + fsize * wspace * (ncols - 1), nrows * fsize)
)
plt.subplots_adjust(wspace=wspace)

# This produces two Axes objects in a single Figure
print("axes:", axs)

blanks = [None] * 3

# change the filtering to adata.obs.condition_status for disease subtypes, or to adata.obs.disease_corrected for simply healthy versus crohns
# use gene_symbols argument to tell scanpy to plot adata.var[gene_symbols] rather than the default adata.var index
alpha_var = 1
s_var = 3

for row_idx, (gene, titles) in enumerate(d_genes.items()):
    print(row_idx)
    print(gene)
    print(titles)
    print("=======")

    adata_subset = adata_subset[adata_subset.obs[gene].argsort()].copy()

    for col_idx, t in enumerate(titles):
        print(col_idx)
        print(t)
        print("------\n")

        ax = axs[row_idx, col_idx]
        # Select the cell types that are not shared in the TenK10K dataset. 
        # If you want to plot ALL, just change adata_subset to adata
        sc.pl.umap(adata_subset[adata_subset.obs.condition_status == t], 
        # sc.pl.umap(adata[adata.obs.condition_status == t], 
                    color=gene, 
                    gene_symbols = "gene_symbols",
                    groups = groups,
                    use_raw = False, # not None!!! this means normalised counts in X are plotted
                    frameon=False,
                    vmax="7",
                    vmin="0",
                    title = t if row_idx == 0 else "", 
                    legend_fontsize=8,
                    legend_fontoutline=0,
                    add_outline=True,        # Adds a nice 'finished' look to clusters
                    outline_width=(0.3, 0.01),
                    alpha = alpha_var,
                    s=s_var,
                    show=False,
                    legend_loc='on data',   # Puts labels over clusters, doesn't seem to work ... 
                    ax = ax)
        # the first Axes object in the row
        if col_idx == 0:
            # We disabled axis drawing in UMAP to have plots without background and border
            # so we need to re-enable axis to plot the ylabel
            ax.axis("on")
            ax.tick_params(
                top="off",
                bottom="off",
                left="off",
                right="off",
                labelleft="off",
                labelbottom="off",
            )
            ax.set_ylabel(gene + "\n", rotation=90, fontsize=16)
            ax.set(frame_on=False)

#description = "mr_deg_interesting_genes_adata_subset_with_labels"
description = "interesting_genes_2_adata_subset_with_labels_ordered"
plt.savefig(f"{figs}/umap_counts_{description}.png", bbox_inches='tight', dpi=300)

