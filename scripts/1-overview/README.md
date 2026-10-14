# Landscape of causal associations between gene expression and complex traits across immune cells

These scripts summarise the scope of the MR analysis: how many genes and traits were tested, how many associations passed local FDR < 0.05 in each cell type and trait, and how the results compare with an earlier OneK1K MR study of coronary artery disease. They produce Figure 1, Supplementary Figures 1–6 and Supplementary Tables 1–3. All scripts run from the repository root.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `study_design.R` | Fig 1b–c panels; Supp Tables 1–2; calls Typst to build Fig 1 | `make fig-m1`, `make tab-s1`, `make tab-s2` | `preprocess_results.R` objects; `metadata/table_column_names.xlsx` | `figures/panels/1-overview/panel_{b,c}.svg`, `figures/fig-m1.{pdf,png}`, `tables/` (Supp Tables 1–2) |
| `fig1.typ` | Fig 1 (assembly of a–c) | `make fig-m1` | `figures/panels/1-overview/infographics_v5.pdf`, `panel_b.svg`, `panel_c.svg` | `figures/fig-m1.{pdf,png}` |
| `mr_heatmap_celltype_diseases.R` | Supp Fig 1 | `make fig-s1` | v5 parquet, MAGMA parquet, `metadata/` | `figures/fig-s1.{pdf,png}` |
| `mr_heatmap_celltype_biotraits.R` | Supp Fig 2 | `make fig-s2` | v5 parquet, MAGMA parquet, `metadata/` | `figures/fig-s2.{pdf,png}` |
| `fig_smr_overview_scatter.R` | Supp Fig 3 (tests vs associations per cell type); Supp Fig 9 (MAGMA vs MR gene counts per trait) | `make fig-s3`, `make fig-s9` | `preprocess_results.R` objects | `figures/fig-s3.{pdf,png}`, `figures/fig-s9.{pdf,png}` |
| `celltype_exclusive_assoc.R` | Supp Fig 4 | `make fig-s4` | `preprocess_results.R` objects | `figures/fig-s4.{pdf,png}` |
| `compare_ray_tenk10k.R` | Supp Fig 5 (CAD MR P values vs Ray et al.); Supp Fig 6 (effect estimates) | `make fig-s5`, `make fig-s6` | `preprocess_results.R` objects; `resources/misc/onek1k_cad_mr_ray2025.txt` | `figures/fig-s5.{pdf,png}`, `figures/fig-s6.{pdf,png}`, `figures/panels/1-overview/compare_ray_tenk10k_overlap.{pdf,png}` (not a manuscript figure) |
| `tab_smr_overview.R` | Supp Table 3 | `make tab-s3` | `preprocess_results.R` objects; `metadata/table_column_names.xlsx` | `tables/` (Supp Table 3) |

## Run order

The scripts are independent of each other. They all source `scripts/0-preprocess/preprocess_results.R` or read `results/preprocessed/tenk10k_phase1.v5.parquet.gz` directly. Figure 1 is compiled by `study_design.R` once it has written the panel SVGs:

```bash
Rscript scripts/1-overview/study_design.R
# or, once the panels exist:
typst compile --root . scripts/1-overview/fig1.typ figures/fig-m1.pdf
typst compile --root . --ppi 300 scripts/1-overview/fig1.typ figures/fig-m1.png
```

## Restricted inputs

None of these scripts read raw single-cell data. Figure 1a needs the BioRender infographic `figures/panels/1-overview/infographics_v5.pdf`, which is not tracked and must be supplied separately. Supplementary Figures 5–6 need the Ray et al. (2025) OneK1K CAD MR results at `resources/misc/onek1k_cad_mr_ray2025.txt`.
