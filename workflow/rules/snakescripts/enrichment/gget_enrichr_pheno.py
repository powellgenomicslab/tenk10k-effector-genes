import gget
import pandas as pd
from pathlib import Path
from gget.compile import PACKAGE_PATH

X = snakemake.wildcards
INPUT = snakemake.input
OUTPUT = snakemake.output
PARAMS = snakemake.params

# interactive:
# X = {'study': 'tenk10k_phase1',
#      'q_thresh': 0.05,
#      'pheno': 'crohns',
#      'heidi_thresh': 0.01}

# INPUT = dict(
#         gene_universe = f"results/enrichment/{X['study']}/gene_universe.txt",
#         msmr_sig = f"results/enrichment_pheno/{X['study']}~q_{X['q_thresh']}~heidi_{X['heidi_thresh']}/gene_celltype/{X['pheno']}.txt",
#         # full gene_set list: https://maayanlab.cloud/Enrichr/#libraries
#         gene_set = "resources/misc/enrichr.gene_set.txt"
# )

# OUTPUT = dict(pheno = f"results/enrichment_pheno/{X['study']}~q_{X['q_thresh']}~heidi_{X['heidi_thresh']}/enrichr/{X['pheno']}.tsv")
# PARAMS = dict(min_gene = 5)

def read_file_to_list(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

df_msmr = pd.read_csv(INPUT['msmr_sig'], sep = " ", names = ["biosample", "gene"])

gene_universe = read_file_to_list(INPUT['gene_universe'])
gene_sets = read_file_to_list(INPUT['gene_set'])

# default gget gene list
with open(f"{PACKAGE_PATH}/constants/enrichr_bkg_genes.txt") as f:
    bg_genes = f.read().splitlines()

biosamples = df_msmr['biosample'].unique().tolist()
dfs_enrich = []

p = X['pheno']
for b in biosamples:
    genes = df_msmr.query("biosample == @b")['gene'].unique().tolist()
    genes = [g for g in genes if g in bg_genes]
    
    # skip if length of genes < PARAMS['min_gene']
    if len(genes) < PARAMS['min_gene']:
        print(f"Skipping enrichment for {p} - {b}: less than {PARAMS['min_gene']} genes")
        continue
    # Enrichr query
    print(f"Enrichment test for {p} - {b} with {len(genes)} genes")
    # run enrichr per gene set, combine into a single dataframe
    df_enrich = pd.concat(
        [
            gget.enrichr(genes, gs, species='human', background_list=bg_genes)
            for gs in gene_sets
        ], ignore_index=True
    )
    df_enrich['phenotype'] = p
    df_enrich['biosample'] = b
    dfs_enrich.append(df_enrich)

# Save results
# exit if dfs_enrich is empty
Path(OUTPUT['pheno']).parent.mkdir(parents=True, exist_ok=True)

if not dfs_enrich:
    print(f"No enrichment results for {p}, writing empty file.")
    Path(OUTPUT['pheno']).touch()
else:
    pd.concat(dfs_enrich, ignore_index=True) \
        .to_csv(OUTPUT['pheno'], sep='\t', index=False)

