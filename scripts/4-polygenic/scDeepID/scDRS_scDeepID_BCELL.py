# Import required libraries
import sys
sys.path.append('/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation')
import scDeepID_minimal2 as scDeepID
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import torch
torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# Load
attn_adata = sc.read("/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/attention_umap.h5ad")
pathway_score = sc.read("/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/pathway_meanpooling.h5ad")
pathway_score_copy = pathway_score.copy()
sc.pp.normalize_total(pathway_score_copy)
sc.pp.scale(pathway_score_copy, max_value=10)
pathway_score_copy.obsm['X_umap'] = attn_adata.obsm['X_umap']

#### SMR integration
# Set up custom font
import matplotlib
import matplotlib.font_manager as font_manager
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
# load smr results, integrate with scdeepid results
df_smr = pd.read_parquet('/g/data/fy54/analysis/tenk10k-causal/results/aggregate/tenk10k_phase1.scdrs.cell_score.tsv.parquet.gz')
df_smr_stat = pd.read_parquet('/g/data/fy54/analysis/tenk10k-causal/results/aggregate/tenk10k_phase1.scdrs.cell_mcp.tsv.parquet.gz')
# filter to B cells
b_cell_types = ['B_naive','B_intermediate', 'B_memory']
df_smr_subset = df_smr[['index' ,'cell_type', 'umap_1', 'umap_2', 'sle']].copy()
df_smr_subset = df_smr_subset[df_smr_subset['cell_type'].isin(b_cell_types)]
df_smr_stat_subst = df_smr_stat[['index', 'cell_type', 'umap_1', 'umap_2', 'sle']].copy()
df_smr_stat_subst = df_smr_stat_subst[df_smr_stat_subst['cell_type'].isin(b_cell_types)]
# create sle labels
common_cells = list(set(pathway_score_copy.obs_names) & set(df_smr_subset['index']))
print(f"Number of common cells: {len(common_cells)}")
pathway_score_subset = pathway_score_copy[common_cells].copy()
df_smr_indexed = df_smr_subset.set_index('index')
df_smr_stat_indexed = df_smr_stat_subst.set_index('index')
pathway_score_subset.obs['sle'] = df_smr_indexed.loc[pathway_score_subset.obs_names, 'sle'].values
pathway_score_subset.obs['sle_pval'] = df_smr_stat_indexed.loc[pathway_score_subset.obs_names, 'sle'].values
print(f"Subsetted shape: {pathway_score_subset.shape}")
print(f"SLE scores added: {'sle' in pathway_score_subset.obs.columns}")
print(f"SLE score range: {pathway_score_subset.obs['sle'].min():.3f} to {pathway_score_subset.obs['sle'].max():.3f}")
pathway_score_subset.obs['umap_1_scdrs'] = df_smr_indexed.loc[pathway_score_subset.obs_names, 'umap_1'].values
pathway_score_subset.obs['umap_2_scdrs'] = df_smr_indexed.loc[pathway_score_subset.obs_names, 'umap_2'].values

### balloon plot
from scipy import stats
# Create SLE status column
pathway_score_subset.obs['sle_status'] = ['SLE' if pval < 0.05 else 'Non-SLE' for pval in pathway_score_subset.obs['sle_pval']]
from statsmodels.stats.multitest import multipletests
# 1. Calculate t-test for all pathways
p_values = []
pathway_names = []
sle_means = []
non_sle_means = []
for pathway in pathway_score_subset.var_names:
    # Extract data
    sle_data = pathway_score_subset[pathway_score_subset.obs['sle_status'] == 'SLE', pathway].X.flatten()
    non_sle_data = pathway_score_subset[pathway_score_subset.obs['sle_status'] == 'Non-SLE', pathway].X.flatten()
    # T-test for SLE > Non-SLE
    t_stat, p_val = stats.ttest_ind(sle_data, non_sle_data, alternative='greater')
    p_values.append(p_val)
    pathway_names.append(pathway)
    sle_means.append(sle_data.mean())
    non_sle_means.append(non_sle_data.mean())
# Calculate FDR
_, fdr_values, _, _ = multipletests(p_values, method='fdr_bh')
# 2. Convert to -log10(p) and rank
neg_log_p = -np.log10(np.array(p_values))
results_df = pd.DataFrame({
    'pathway': pathway_names,
    'neg_log_p': neg_log_p,
    'p_value': p_values,
    'fdr': fdr_values,
    'sle_mean': sle_means,
    'non_sle_mean': non_sle_means
})
results_df_raw = results_df.copy()
results_df_raw = results_df_raw.sort_values('neg_log_p', ascending=False)
# Filter: FDR < 0.05 AND at least one group should have mean > 0.2
results_df = results_df[(results_df['fdr'] < 0.05) & ((results_df['sle_mean'] > 0.2) | (results_df['non_sle_mean'] > 0.2))]
# Sort by -log10(p)
results_df = results_df.sort_values('neg_log_p', ascending=False)
# Print FDR summary
print(f"Total pathways tested: {len(p_values)}")
print(f"FDR < 0.05: {sum(fdr_values < 0.05)}")
print(f"FDR < 0.05 & mean > 0.2: {len(results_df)}")
# 3. Create dumbbell plot for top FDR-significant pathways
n_pathways = min(50, len(results_df)) 
top_pathways = results_df.head(n_pathways)
fig, ax = plt.subplots(figsize=(10, 8)) 
for i, (idx, row) in enumerate(top_pathways.iterrows()):
    # Reverse the y-position
    y_pos = n_pathways - 1 - i
    ax.plot([row['non_sle_mean'], row['sle_mean']], [y_pos, y_pos], 'k-', alpha=0.5)
    ax.scatter(row['non_sle_mean'], y_pos, color='lightblue', s=60, label='Non-SLE' if i==0 else "")
    ax.scatter(row['sle_mean'], y_pos, color='lightcoral', s=60, label='SLE' if i==0 else "")
ax.set_yticks(range(n_pathways))
reversed_labels = top_pathways['pathway'].str.replace('GOBP_', '').str.replace('_', ' ').tolist()[::-1]
ax.set_yticklabels(reversed_labels, fontsize=8)
ax.set_xlabel('Mean Pathway Score')
ax.set_title(f'Top {n_pathways} FDR-Significant Pathways (FDR < 0.05, mean > 0.2)', fontsize=14)
ax.legend()
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig('/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/smr_plot/top50_fdr_pathways.pdf', dpi=300,bbox_inches='tight')  
plt.show()
print(f"\nTop {n_pathways} FDR-significant pathways:")
for i, row in top_pathways.iterrows():
    print(f"{row['pathway']}: p={row['p_value']:.3e}, FDR={row['fdr']:.3e}, -log10(p)={row['neg_log_p']:.2f}")


# Plot pathway umap
pathways_to_plot = top_pathways[0:5]['pathway'].tolist()  
pathways_to_plot_umap = pathways_to_plot.copy()
pathways_to_plot_umap.insert(0, 'celltype')
fig = sc.pl.umap(pathway_score_subset, color=pathways_to_plot_umap,
                   cmap='magma_r', vmax=5, size=3, ncols=5,
                   return_fig=True, show=False)
for ax in fig.get_axes():
      for artist in ax.collections:
          artist.set_rasterized(True)
fig.savefig('/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/smr_plot/pathway_umap_custom.pdf', dpi=300, bbox_inches='tight')
plt.close(fig)
for pathway in pathways_to_plot_umap:
      if pathway != 'celltype':
          max_val = pathway_score_subset[:, pathway].X.max()
          print(f"{pathway}: {max_val:.3f}")

# Plot violin + boxplot
fig, axes = plt.subplots(1, 5, figsize=(20, 4))
colors = ['lightblue', 'lightcoral']  # Non-SLE: blue, SLE: red
for i, pathway in enumerate(pathways_to_plot):
    # Extract data
    sle_data = pathway_score_subset[pathway_score_subset.obs['sle_status'] == 'SLE', pathway].X.flatten()
    non_sle_data = pathway_score_subset[pathway_score_subset.obs['sle_status'] == 'Non-SLE', pathway].X.flatten()
    # Calculate means
    sle_mean = sle_data.mean()
    non_sle_mean = non_sle_data.mean()
    # Perform one-sided t-test (SLE > Non-SLE)
    t_stat, p_value_greater = stats.ttest_ind(sle_data, non_sle_data, alternative='greater')
    # Create violin plot
    parts = axes[i].violinplot([non_sle_data, sle_data], positions=[1, 2], showmeans=False, showmedians=False, showextrema=False)
    for j, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[j])
        pc.set_alpha(0.6)
    # Overlay boxplot without outliers
    bp = axes[i].boxplot([non_sle_data, sle_data], positions=[1, 2],
                        widths=0.15, showfliers=False,
                        patch_artist=True,
                        boxprops=dict(facecolor='white', color='black'),
                        medianprops=dict(color='black', linewidth=2))
    axes[i].set_xticks([1, 2])
    axes[i].set_xticklabels(['Non-SLE', 'SLE'])
    axes[i].set_title(pathway.replace('GOBP_', '').replace('_', ' '), fontsize=8)
    axes[i].set_ylabel('Pathway Score')
    direction = "↑" if sle_mean > non_sle_mean else "↓"
    axes[i].text(0.5, 0.95, f'SLE {direction}\np = {p_value_greater:.3e}',
                transform=axes[i].transAxes,
                ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.tight_layout()
plt.savefig('/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/smr_plot/top5violin.pdf', dpi=300,bbox_inches='tight')
plt.show()

results_df_raw.to_csv("/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/smr_plot/pathway_test.csv", index=False)
top5_pathways = results_df.head(5)['pathway'].tolist()
plt.figure(figsize=(10, 4))  
colors = ['lightcoral', 'lightskyblue', 'lightgreen', 'peachpuff', 'plum']
for pathway, color in zip(top5_pathways, colors):
    sle_scores = pathway_score_subset.obs['sle'].values
    pathway_scores = pathway_score_subset[:, pathway].X.flatten()
    # Calculate 10 quantiles
    quantiles = pd.qcut(sle_scores, q=10, labels=['Q1', 'Q2', 'Q3', 'Q4', 'Q5',
                                                   'Q6', 'Q7', 'Q8', 'Q9', 'Q10'])
    # Calculate means per quantile
    means = []
    for q in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'Q9', 'Q10']: means.append(pathway_scores[quantiles == q].mean())
    plt.plot(range(1, 11), means, marker='o', markersize=8, color=color, linewidth=2, label=pathway.replace('GOBP_', '').replace('_', ' '))
plt.xlabel('scDRS Score Quantiles', fontsize=14)
plt.ylabel('Mean Pathway Score', fontsize=14)
plt.title('Top 5 Pathways by SLE Score (10 Quantiles)', fontsize=16)
plt.xticks(range(1, 11), ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'Q9', 'Q10'])
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('/g/data/ei56/projects/lawhua/proj_scDeepID/minimal_ablation/exp_meanpooling/tk_bcell_final/smr_plot/top5_quintile.pdf', dpi=300, bbox_inches='tight')
plt.show()