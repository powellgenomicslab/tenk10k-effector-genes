import scanpy as sc
import pandas as pd
from pathlib import Path
import yaml

INPUT = snakemake.input
OUTPUT = snakemake.output


# Load data
adata = sc.read(INPUT['h5ad'])

# Set index to hgnc symbol
adata.var.set_index(['gene_name'], inplace=True, drop=False)

# replace to get the raw count
adata = adata.raw.to_adata()

# rewrite adata object
adata.write_h5ad(OUTPUT['prep_h5ad'])