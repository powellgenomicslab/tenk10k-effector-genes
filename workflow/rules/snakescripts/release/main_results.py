"""Write the main results of the study as a standalone, self-describing file.

Input: results/preprocessed/tenk10k_phase1.v5.parquet.gz, the integrated MR,
sensitivity and colocalisation table that the figure and table scripts read.
The release copy keeps every column, renamed and ordered as in
workflow/config/main_results_columns.tsv (the names follow the Hugging Face
release where it has the same column). Each column's description is stored in
the Parquet field metadata and the dataset description in the schema metadata;
the same table is written out as a Markdown column dictionary. Rows are sorted
by trait, cell type and gene.

Run by rule main_results_release in workflow/rules/misc.smk.
"""

import json

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

INPUT = snakemake.input.results
OUT_PARQUET = snakemake.output.parquet
OUT_MD = snakemake.output.md
COLUMNS = snakemake.params.columns
META = dict(snakemake.params.meta)

cols = pd.read_csv(COLUMNS, sep="\t")
df = pd.read_parquet(INPUT)

missing = set(df.columns) - set(cols["source"])
extra = set(cols["source"]) - set(df.columns)
if missing or extra:
    raise ValueError(f"column table does not match the input: missing {missing}, extra {extra}")

df = df[list(cols["source"])].rename(columns=dict(zip(cols["source"], cols["column"])))
for c in ("cell_type", "major_cell_type", "pheno_cat"):
    df[c] = df[c].astype(str)
df = df.sort_values(["pheno_id", "cell_type", "ensembl_gene_id"], kind="stable").reset_index(drop=True)

table = pa.Table.from_pandas(df, preserve_index=False)
described = dict(zip(cols["column"], cols["description"]))
schema = pa.schema(
    [f.with_metadata({"description": described[f.name]}) for f in table.schema],
    metadata={
        **{k: str(v) for k, v in META.items()},
        "n_rows": str(len(df)),
        "columns": json.dumps(described),
    },
)
pq.write_table(table.cast(schema), OUT_PARQUET, compression="zstd")

with open(OUT_MD, "w") as fh:
    fh.write(f"# {META['title']}\n\n{META['description']}\n\n")
    fh.write(f"File: `{OUT_PARQUET.split('/')[-1]}`; {len(df):,} rows (one per cell type, trait and gene "
             f"tested) and {df.shape[1]} columns. Version: {META['version']}. Licence: {META['licence']}. "
             f"Data DOI: {META['doi']}.\n\n")
    fh.write("Each column's description is also stored in the Parquet file's field metadata, "
             "and this dataset description in its schema metadata.\n")
    for group, rows in cols.groupby("group", sort=False):
        fh.write(f"\n## {group}\n\n| Column | Description |\n|---|---|\n")
        for _, r in rows.iterrows():
            fh.write(f"| `{r['column']}` | {r['description']} |\n")

print(f"wrote {OUT_PARQUET}: {len(df):,} rows, {df.shape[1]} columns")
