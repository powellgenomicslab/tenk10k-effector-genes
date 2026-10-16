# Releasing the code and the data

The code and the data are archived as two Zenodo records, each with its own DOI:

- **Code**: this repository. Each GitHub release is archived by Zenodo's GitHub
  integration, using the metadata in `.zenodo.json`.
- **Data**: the bundle of figure and table inputs (`make bundle`), uploaded as
  a Zenodo dataset with the metadata in `zenodo/data.json`.

They are linked both ways: `make data` downloads the bundle from the data
record and checks every file against its SHA-256, and each record's metadata
cites the other's DOI. A reader who clones a release and runs `make data` gets
exactly the inputs that release was tested with.

## Once

1. The `creators` in `.zenodo.json` (code) and `zenodo/data.json` (data)
   follow the manuscript's author list and affiliations. Add ORCIDs
   (`"orcid"` per creator) if wanted. The data DOI reserved for the first
   release is 10.5281/zenodo.23406585 (record 23406585).
2. A repository owner enables the repository in Zenodo under
   *Account → GitHub*.

## Each release

1. **Reserve the data DOI.** Start a new upload on Zenodo (or a new version of
   the existing data record) and click *Reserve DOI*. Note the DOI and the
   record number.
2. **Point the code at the data.** In `config/figures.yaml`, set
   ```yaml
   bundle_url: "https://zenodo.org/records/<record>/files/tenk10k_figure_inputs.tar.gz?download=1"
   ```
   Put the data DOI in the README's *Data Availability* section and in the
   `related_identifiers` of `.zenodo.json` (relation `isSupplementedBy`).
   Commit.
3. **Build the bundle** on a compute node, from the commit to be released:
   ```bash
   make link-data DATA_DIR=/g/data/fy54/analysis/tenk10k-causal
   make bundle BUNDLE=tenk10k_figure_inputs.tar.gz
   ```
   This writes `INVENTORY.tsv` and stops if any file lacks an entry in
   `config/figure_inputs_annotation.tsv` or if a checksum does not match.
4. **Test the bundle as a reader.** In a fresh clone of the commit, temporarily
   set `bundle_url` to `file:///path/to/tenk10k_figure_inputs.tar.gz`, then run
   `make data`, `make dry-run` (no missing inputs), and build a few items on a
   compute node.
5. **Release the code.** Tag and publish a GitHub release (e.g. `v4.1.0`).
   Zenodo archives it and mints the code DOI.
6. **Publish the data.** Upload the main results,
   `tenk10k_phase1_effector_genes.v5.parquet` and its column dictionary
   `tenk10k_phase1_effector_genes.v5.md` (rule `main_results_release`, written
   to `results/release/` in the working repository); then
   `tenk10k_figure_inputs.tar.gz`, and also `INVENTORY.tsv` and
   `scripts/0-preprocess/BUNDLE_README.md` (as `README.md`) so the contents can
   be read without downloading the archive.
   Complete the metadata from `zenodo/data.json`, add the code DOI and the
   release tag under `related_identifiers` (relation `isSupplementTo`), and
   publish.
7. **Check the link.** In a fresh clone of the release tag, run `make data` and
   build one item.

## Licences

The code is under the MIT licence (`LICENSE`). The data record is under
CC BY-SA 4.0; files from third parties keep their own terms, which
`INVENTORY.tsv` lists file by file.
