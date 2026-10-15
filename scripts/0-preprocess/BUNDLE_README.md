# TenK10K Phase 1 effector genes: inputs for the figures and tables

This archive holds the summary-level data needed to rebuild every figure and
supplementary table of *Single-cell genetics identifies cell-type-specific
effector genes across complex traits and diseases* (Henry, Senabouth,
Tyebally et al.; medRxiv doi:10.1101/2025.08.28.25334614). The code is at
<https://github.com/powellgenomicslab/tenk10k-effector-genes>.

## Using the archive

The repository's `make data` downloads this archive, checks every file against
`figure_inputs.sha256` and unpacks it in place. To do this by hand, unpack the
archive in the root of a clone of the repository:

```bash
git clone https://github.com/powellgenomicslab/tenk10k-effector-genes.git
cd tenk10k-effector-genes
tar -xzf /path/to/tenk10k_figure_inputs.tar.gz
shasum -a 256 -c figure_inputs.sha256
make fig-m2
```

The repository's README describes the software and the full set of `make`
targets.

## Contents

| File | Contents |
|---|---|
| `INVENTORY.tsv` | One row per data file: path, size, SHA-256, category, description, the display items that read it, source and licence |
| `figure_inputs.sha256` | SHA-256 checksums of the data files |
| `results/` | Results of the analysis pipelines: Mendelian randomisation, colocalisation, MAGMA, scDRS, enrichment and Open Targets queries, and the locus data behind the regional plots |
| `resources/` | Reference tables and extracts of published tables used by the Crohn's disease case study and the enrichment analyses |

`INVENTORY.tsv` lists the categories:

- **TenK10K results**: produced by this study.
- **Summary-statistics excerpts**: published GWAS and eQTL summary statistics
  for the few loci that are plotted.
- **Public-data queries and public references**: Open Targets, ChEMBL,
  STRING, MSigDB and Gene Ontology tables.
- **Third-party tables and extracts**: published supplementary tables. Two
  of them, Kong et al. 2023 and Liu et al. 2023, are included only as the
  extracts the analysis uses.

## Not included

- Individual-level single-cell expression data and genotypes of the TenK10K
  donors. The LD values in `*.ld` files and in the locus objects are r²
  statistics computed from the genotypes; no genotypes are included. In the
  cell-level scDRS files, each cell ID is the 10x barcode followed by an
  anonymous code in place of the sequencing pool.
- Full GWAS and eQTL summary statistics. These are available from their
  original sources (see `INVENTORY.tsv`). The TenK10K sc-eQTL summary
  statistics are released with Cuomo et al. 2025
  (doi:10.1101/2025.03.20.25324352).
- The single-cell data of Kong et al. 2023 and of the AIFI Immune Health
  Atlas (Gong et al. 2025), and the original supplementary tables of Kong et
  al. 2023 and Liu et al. 2023.

## Licences

Each file's licence is given in the `licence` column of `INVENTORY.tsv`.
The archive is released under CC BY-SA 4.0, which also covers the files
produced by this study. Files derived from ChEMBL are under CC BY-SA 3.0.
Open Targets Platform outputs are under CC0. Third-party tables and
summary-statistics excerpts keep the terms of their sources, and should be
cited as given in the `source` column. The STRING enrichment results include
KEGG pathway terms (rows with source `KEGG`), which remain under KEGG's terms
(<https://www.kegg.jp/kegg/legal.html>).

## Citation

Henry A, Senabouth A, Tyebally R, et al. Single-cell genetics identifies
cell type-specific causal mechanisms in complex traits and diseases.
medRxiv 2025.08.28.25334614 (2025). doi:10.1101/2025.08.28.25334614
