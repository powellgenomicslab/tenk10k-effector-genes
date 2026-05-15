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
# Plot Genes of Interest
# Either by Type (adata.obs.condition_status, nCat = 3) or Disease (adata.obs.disease_corrected, nCat = 2). See prepare_data_for_MR_comparison script.
################################################################################################################################################################################################################################

diseases = ['Healthy', 'Crohns Disease']
types = ['Healthy Colon', 'Non-inflamed Crohns Disease Colon', 'Inflamed Crohns Disease Colon']

# Set the axes of the final figure = columns are the conditions, rows are the genes. For diseases ncols = 2, for types, ncol =3
# adjust wspace for space between figures 

ncols = 3
nrows = 2
fsize = 4
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
# d_genes = {'AIF1': types, #MATCHED DIR
#             'PTGER4': types, #MATCHED DIR
#             'CDC42SE2': types, #MATCHED DIR
#             'GPX1': types, # MIXED MATCHED DIR
#             'TNFRSF18': types, # False DIR
#             'BIN1': types} # False DIR

d_genes = {'NCF4': diseases,
            'ZBTB38': diseases}
d_genes = {'NCF4': types,
            'ZBTB38': types}

# change the filtering to adata.obs.condition_status for disease subtypes, or to adata.obs.disease_corrected for simply healthy versus crohns
# use gene_symbols argument to tell scanpy to plot adata.var[gene_symbols] rather than the default adata.var index

for row_idx, (gene, titles) in enumerate(d_genes.items()):
    print(row_idx)
    print(gene)
    print(titles)
    print("=======")

    for col_idx, t in enumerate(titles):
        print(col_idx)
        print(t)
        print("------\n")

        ax = axs[row_idx, col_idx]
        #sc.pl.umap(adata[(adata.obs.disease_corrected == t) & (~adata.obs['major_cell_type'].isin(['Immune Cycling cells', 'Macrophage', 'Mast', 'Intraepithelial NK-like']))], 
        sc.pl.umap(adata[adata.obs.condition_status == t], 
                    color=gene, 
                    gene_symbols = "gene_symbols",
                    use_raw = None,
                    frameon=False,
                    vmax="6",
                    vmin="0",
                    title = t if row_idx == 0 else "", 
                    legend_fontsize=6,
                    legend_fontoutline=0,
                    s=20,
                    show=False,
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
    
plt.savefig(f"{figs}/6_GENES_all_cell_types_disease_subtypes.png", bbox_inches='tight', dpi=300)

