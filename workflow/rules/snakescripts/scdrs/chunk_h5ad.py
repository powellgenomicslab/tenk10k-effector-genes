import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

STUDY = "immune_atlas"
PARAMS = {"n_chunks": 10}
H5AD = f"resources/scdrs/h5ad/{STUDY}.h5ad"

adata = sc.read_h5ad(H5AD)

adata.chunk_X(select = 10000)