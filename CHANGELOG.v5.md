# `tenk10k_phase1.v5.parquet.gz` — changelog

Preprocessed integration of SMR/MSMR, IVW-LD, MRlink2, coloc and mvcoloc evidence
for TenK10K Phase 1. Produced by the `preprocess_results` rule
(`workflow/rules/misc.smk`) from
`workflow/rules/snakescripts/preprocess/tenk10k_phase1.v5.R`.

Column definitions follow `results/preprocessed/tenk10k_phase1.v2.md`; the schema
is unchanged from v4 apart from `n_snps_smr` (below).

---

## Input sources

| evidence | source |
|---|---|
| MSMR | `results/sensitivity/smr/tenk10k_phase1/tenk10k_phase1_sensitivity.msmr.parquet.gz` |
| IVW-LD | `results/sensitivity/ivw-ld/robust/tenk10k_phase1_sensitivity.ivw-ld.v2.parquet.gz` |
| MRlink2 | `results/sensitivity/mrlink2/tenk10k_phase1_sensitivity.mrlink2.v2.parquet.gz` |
| coloc | `results/aggregate/coloc/tenk10k_phase1.coloc.v3.parquet.gz` |
| mvcoloc | `results/aggregate/coloc/tenk10k_phase1.mvcoloc.parquet.gz` (unchanged) |

---

## Changes relative to v4

### 1. New column `n_snps_smr`

Number of SNPs actually used as instruments in multi-SNP SMR per
`(biosample, phenotype, gene)`, counted from
`results/smr/tenk10k_phase1/*/<trait>/all_chr.snps4msmr.list`. Lets downstream
filtering distinguish single- from multi-instrument tests.

### 2. IVW-LD and MRlink2 moved to the `v2` results

v4 used the original `…ivw-ld.parquet.gz` / `…mrlink2.parquet.gz`; v5 uses the
`v2` outputs. The IVW source is specifically the **`robust`** variant. Paths are
now written out in full rather than via the `results/aggregate/` symlinks, so the
source is explicit in the script.

### 3. MSMR results for `asd` replaced from per-biosample SMR output

`asd` SMR was rerun after the sensitivity MSMR parquet was built, so v5 reads
`results/smr/tenk10k_phase1/*/asd/all_chr.msmr` directly and substitutes those
rows (pre-existing v5 behaviour, retained).

### 4. coloc moved from v1 to **v3** — the substantive change

v4 consumed `tenk10k_phase1.coloc.parquet.gz` (v1). v5 consumes **v3**, which
fixes two independent classes of defect. Row counts:

| artifact | rows | note |
|---|---|---|
| v1 (used by v4) | 11,594,987 | had a junk `V13` column; 894 rows with `NA` SNP count |
| v2 | 11,930,730 | aggregation fixed only; same underlying coloc results |
| **v3 (used by v5)** | **15,434,261** | + coloc recomputed for 34 traits |

**4a. Aggregation fixes (v1 → v2), +335,743 rows.** The previous reader used
`fread(fill = TRUE)`, which silently stops at the first row wider than the header
and discards the remainder — `ASDC/alzheimers` returned 317 of its 653 rows with
only a warning. ~80 aggregated files were affected, a consequence of
`concat_coloc_chr` having emitted only the first chromosome's header while
appending rows of any width. The SNP count could also land in an unnamed 13th
field that the old coalesce could not see (the source of v1's 894 `NA`s and its
`V13` column). The reader now parses line-by-line, cannot truncate, coalesces the
count from all three possible columns, and errors if any `V\d+` column survives or
any successful row lacks a count.

**4b. coloc recomputed for 34 traits (v2 → v3), +3,503,531 rows.** Two GWAS
data-quality defects made `coloc.abf` fail for *every* gene whose ±100 kb window
contained one bad SNP, so affected traits lost whole chromosomes:

- `b = 0` with `se = 0` (unreported effects) → `z = 0/0 = NaN`, so `min(p)` is `NA`
  and `check_dataset()` aborts with "missing value where TRUE/FALSE needed".
  `sle` chr1: 48,700 of 581,834 SNPs.
- duplicate SNP IDs — multi-allelic sites collapsed onto one ID — which multiply
  the eQTL rows through the keyed join, so `check_dataset()` rejects dataset 1
  with "duplicated snps found". `pancreatic` chr3: 25,833.

`run_coloc_test.R` now drops non-finite/zero-`se` SNPs and ambiguous duplicate
IDs before the join, logging what it removed.

Recovery, as % of `asthma` (a clean reference trait):

| pheno | v2 | v3 |
|---|---|---|
| `sle`, `cmelanoma`, `pancreatic` | 0.5% | 100% |
| `breast_cimba` | 1.1% | 100% |
| `ra` | 3.0% | 100% |
| `asd` | 8.4% | 100% |
| `parkinsons` | 10.6% | 100% |
| `psoriasis` | 11.8% | 100% |
| `SkBmd` | 24.0% | 100% |
| `prostate` | 53.7% | 100% |
| `t1dm` | 76.7% | 100% |

Lost (biosample × pheno × chromosome) slots fell from **2,326 to 2**; both
remaining belong to `cardio_meta`, which is excluded from
`resources/misc/target_phenotypes.txt`. **No trait now sits below 90% of
`asthma`.** All recomputed traits converge on exactly 154,932 rows — expected,
since the eGene set is phenotype-independent, so once no genes are silently
dropped every trait must yield the same number of gene × biosample tests.

**4c. New `status` column in the coloc source.** `df_blank()` previously returned
*zero* rows on error, so failed genes vanished and were indistinguishable from
genes never attempted. Failures now emit one row carrying the `coloc.abf` error.
In v3, 14,957 rows (0.098%) are labelled failures — `0 (non-NA) cases`, i.e. genes
with no overlapping SNPs. These are the only rows with a `NA` SNP count, and were
previously invisible rather than absent.

### Net effect on the released table

Row count is unchanged at **7,546,764** (the MSMR backbone determines it), but coloc
evidence is now attached to almost all of it:

| | v4 | v5 |
|---|---|---|
| rows with `coloc_pph4` | 5,621,775 (74.5%) | **7,534,793 (99.8%)** |
| traits improved | — | 38 |
| traits unchanged | — | 62 |
| traits worse | — | **0** |

Per-trait coloc annotation rate:

| pheno | v4 | v5 |
|---|---|---|
| `sle` | 0.4% | 100% |
| `cmelanoma` / `pancreatic` | 0.5% | 100% |
| `ra` | 2.8% | 100% |
| `asd` | 7.4% | 100% |
| `asthma` | 30.4% | 100% |
| `prostate` | 52.4% | 100% |
| `t1dm` | 75.6% | 100% |

`asthma` is worth noting because its coloc results were never *computed* wrongly —
it was one of the traits missing ~106,953 rows purely through the aggregation
truncation (§4a). Its jump from 30.4% to 100% therefore comes from the reader fix
alone, not from any recomputation.

### 5. Evidence-flag NA handling fixed

The OR'd evidence criteria produced `NA` instead of `FALSE` whenever one source was
negative and another missing, because R evaluates `FALSE | NA` as `NA`. (`TRUE | NA`
is `TRUE`, so a single positive source always won — that part was correct.)

Since `mvcoloc_pph4` is `NA` for 93.6% of rows — mvcoloc legitimately finds nothing
for most gene × trait pairs — **93% of `coloc` values were `NA`**, of which
7,005,035 had a known `coloc_pph4 < 0.8` and should have been `FALSE`. The
practical consequence: `filter(coloc == FALSE)` returned **420,612** rows instead
of **7,425,647**, because `NA == FALSE` is `NA`. Any downstream count of genes
*without* coloc support was wrong by an order of magnitude.

Each comparison is now wrapped in `%in% TRUE` (maps `NA` to `FALSE`, i.e. a missing
source is non-supporting), with `NA` preserved only where *no* source was assessed,
so "not supported" stays distinguishable from "not tested":

| flag | TRUE before → after | NA before → after |
|---|---|---|
| `coloc` | 109,161 → 109,161 | 7,016,991 → **11,487** |
| `sensitivity` | 6,425,020 → 6,425,020 | 1,006,740 → **107,231** |

No positive call changed, and `max_evidence` is unaffected (`case_when` already
treated an `NA` condition as non-matching, so every MR-positive gene still received
its correct label — `NA` count 7,144,416 before and after).

---

## Caveats

- **v3 coloc is a hybrid.** 34 traits were recomputed; the other 66 carry forward
  prior results. All 66 sit at ≥94% of `asthma` (most ≥98.9%), so they were not
  materially affected, but they have not been recomputed. Rows from those traits
  have `status = NA`.
- **Upstream GWAS defects are not fixed.** `run_coloc_test.R` cleans defensively
  at read time, which protects coloc, but the `.ma` files themselves still carry
  the zero-SE and duplicate-ID SNPs. The MR inputs to this release therefore still
  read them: ~1–2.4% of instruments for `ra`/`cmelanoma`/`pancreatic` are
  arbitrary picks between conflicting duplicate records (some sign-flipped and
  nominally significant), and ~9.3% of `sle` instruments are zero-SE, yielding
  `p_SMR = NA`. Outstanding upstream fixes: exclude multi-allelic sites and join on
  alleles in `format_gwas/ra.R`; add `sort -u` to `mk_tabix_region`; investigate the
  ~706k-SNP zero-SE source family and `ms`.
- **mvcoloc is unchanged** and was not assessed for the same defects.

---

## Supporting changes in this release

- `workflow/rules/coloc.smk` — `concat_coloc_chr` is header-aware and fails on
  mismatch (was `awk 'NR == 1 || FNR > 1'`, which interleaved incompatible
  schemas); `run_coloc` reduced from `threads/ncpus: 4` to `1`, the workload being
  strictly serial (measured 16 CPU-seconds against 39 wall-seconds on 4 cores).
- `workflow/rules/snakescripts/aggregate/coloc_schema.R` — new shared reader.
- `workflow/rules/snakescripts/coloc/run_coloc_test.R` — GWAS cleaning; failure
  rows with `status`; working `log:` capture (`on.exit` under `source()` fires
  immediately, so nothing had ever been written and `logs/coloc/` did not exist);
  and `include == TRUE` when selecting trait metadata, matching the sibling
  scripts — without it the duplicate BMD study rows made `pheno_type` a length-2
  vector and aborted the script under R ≥ 4.2.
- `scripts/diagnostic/` — census and defect-scan tooling used to quantify the
  above (`scan_coloc_results.sh`, `scan_ma_defects*.sh`,
  `classify_missing_slots.R`, `check_mr_instruments.sh`).
- `scripts/aggregate/rebuild_coloc_parquet.{R,sh}` — standalone parquet rebuild
  that does not evaluate the snakemake DAG (which would delete surviving `temp()`
  per-chromosome files).
