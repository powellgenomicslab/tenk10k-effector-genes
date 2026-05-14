import gget
import pandas as pd
from gget.compile import PACKAGE_PATH

X = snakemake.wildcards
INPUT = snakemake.input
OUTPUT = snakemake.output

# interactive:
# X = {'study': 'tenk10k_phase1',
#      'q_thresh': 0.05,
#      'heidi_thresh': 0.01}

# INPUT = dict(
#         gene_universe = f"results/enrichment/{X['study']}/gene_universe.txt",
#         msmr_sig = f"results/aggregate/msmr_sig/{X['study']}~q_{X['q_thresh']}~heidi_{X['heidi_thresh']}.msmr_sig.tsv",
#         # full gene_set list: https://maayanlab.cloud/Enrichr/#libraries
#         gene_set = "resources/misc/enrichr.gene_set.txt"
# )

def read_file_to_list(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

df_msmr = pd.read_csv(INPUT['msmr_sig'], sep = "\t")

gene_universe = read_file_to_list(INPUT['gene_universe'])
gene_sets = read_file_to_list(INPUT['gene_set'])

# default gget gene list
with open(f"{PACKAGE_PATH}/constants/enrichr_bkg_genes.txt") as f:
    bg_genes = f.read().splitlines()


for i in ['biosample', 'phenotype']:
    dfs_enrich = []
    for j in df_msmr[i].unique():
        genes = df_msmr[df_msmr[i] == j]['Gene'].unique().tolist()
        genes = [g for g in genes if g in bg_genes]
        
        # Enrichr query
        print(f"Enrichment test for {j} with {len(genes)} genes")
        # run enrichr per gene set, combine into a single dataframe
        df_enrich = pd.concat(
            [
                gget.enrichr(genes, gs, species='human', background_list=bg_genes)
                for gs in gene_sets
            ], ignore_index=True
        )
        df_enrich[i] = j
        dfs_enrich.append(df_enrich)
    
    # Save results
    pd.concat(dfs_enrich, ignore_index=True) \
        .to_csv(OUTPUT[i], sep='\t', index=False)

