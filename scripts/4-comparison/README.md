# Single-cell eQTL enhances identification of effector genes for complex traits

These scripts compare the single-cell eQTL MR results with two established alternatives: MAGMA gene-based association, and MR using eQTLGen whole-blood instruments. They tally the gene–trait pairs that overlap and those found by only one method, by trait and by the number of cell types. They also draw the two example loci in Fig. 3: a multi-gene eczema locus on chromosome 12, and *BANK1* in SLE.

| Script | Manuscript item(s) | Make target(s) | Key inputs | Outputs |
|---|---|---|---|---|
| `mr_example.R` | Fig. 3 | `make fig-m3` | `scripts/0-preprocess/preprocess_results.R` inputs; `metadata/gencode.v44.gene_type.tsv.gz`; `results/sensitivity/smr/tenk10k_phase1/CD14_Mono/eczema/all_chr.snps4msmr.list`; `results/smr_locus/{tenk10k_phase1,eqtlgen2020}/…`; `results/gwas_locus/{eczema,sle}/…` | `figures/fig-m3.{png,pdf}` |
| `mr_overlap.R` | Extended Data Fig. 4a–b | `make fig-e4` | `scripts/0-preprocess/preprocess_results.R` inputs | `figures/fig-e4.{png,pdf}` |
| `fig_mr_gene_overlap.R` | Supplementary Figs 10 and 11 | `make fig-s10`, `make fig-s11` | `scripts/0-preprocess/preprocess_results.R` inputs | `figures/fig-s10.{png,pdf}`, `figures/fig-s11.{png,pdf}` |
| `tab_mr_overlap.R` | Supplementary Tables 6 and 7 | `make tab-s6`, `make tab-s7` | `scripts/0-preprocess/preprocess_results.R` inputs; `metadata/table_column_names.xlsx` | `tables/tab-s6.tsv`, `tables/tab-s7.tsv` |

Supplementary Fig. 9, which correlates MR and MAGMA gene counts per trait, is drawn by `scripts/1-overview/fig_smr_overview_scatter.R`. That script lives with the Results section 1 scripts because it also draws Supplementary Fig. 3. The script that produced Fig. 3 in the source repository also drew Supplementary Fig. 4; that plot is now in `scripts/1-overview/celltype_exclusive_assoc.R`.

The scripts are independent of one another and can run in any order.

## Restricted inputs

`mr_example.R` reads per-locus eQTL, GWAS and LD dumps under `results/smr_locus/` and `results/gwas_locus/`. The LD matrices come from individual-level TenK10K genotypes and are not part of the public data release, so a third party cannot regenerate the Fig. 3 locus panels. The summary comparisons can be regenerated.
