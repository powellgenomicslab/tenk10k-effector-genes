# Releasing the code and the data

The code and the data are archived as two Zenodo records, each with its own DOI:

- **Code**: this repository. Each GitHub release is archived by Zenodo's GitHub
  integration, using the metadata in `.zenodo.json`. Concept DOI (all versions)
  10.5281/zenodo.23487289; v4.0 is 10.5281/zenodo.23487290.
- **Data**: the bundle of figure and table inputs (`make bundle`), uploaded as
  a Zenodo dataset with the metadata in `zenodo/data.json`. DOI
  10.5281/zenodo.23406585 (files restricted until the manuscript is published).

Both records are in the `powellgenomicslab` community. Their creators, with
every author's affiliations, come from `zenodo/creators.json`, which follows
`resources/authors.docx`.

They are linked both ways: `make data` downloads the bundle from the data
record and checks every file against its SHA-256, and each record's metadata
cites the other's DOI. A reader who clones a release and runs `make data` gets
exactly the inputs that release was tested with.

## Once

1. The `creators` in `.zenodo.json` (code) and `zenodo/data.json` (data)
   follow the manuscript's author list, with each author's first affiliation
   only: both files use Zenodo's legacy schema, which allows one. The full
   lists are in `zenodo/creators.json`. The data DOI reserved for the first
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
5. **Publish the data.** The record holds five files, all in `data/zenodo/`
   once the bundle is built: the main results,
   `tenk10k_phase1_effector_genes.v5.parquet`, and its column dictionary
   `tenk10k_phase1_effector_genes.v5.md` (rule `main_results_release`, written
   to `results/release/` in the working repository); the bundle,
   `tenk10k_figure_inputs.tar.gz`; and `INVENTORY.tsv` and
   `scripts/0-preprocess/BUNDLE_README.md` (as `README.md`), so the contents can
   be read without downloading the archive. Keep the file names: `make data`
   downloads the bundle by name. Upload them with `zenodo/upload_data.py`
   (below) or in the browser, complete the metadata from `zenodo/data.json`,
   set access to public, and submit the draft to the `powellgenomicslab`
   community; it is published when a community curator accepts it. Files cannot
   be changed once published; a change needs a new version, with a new DOI.
6. **Check the link.** In a fresh clone, run `make data` (it downloads from the
   published record and checks every file) and build one item.
7. **Release the code** (can come later). Tag and publish a GitHub release
   from `main` (e.g. `v4.1.0`); Zenodo archives it under the code concept DOI.
   Then restore every author's affiliations on the new version, which
   `.zenodo.json` cannot hold:
   ```bash
   python3 zenodo/set_creators.py --record <new code record ID> --publish
   ```
   The data record already lists the code concept DOI under related works
   (relation "is supplement to"), so it needs no change; editing metadata
   would not create a new version in any case.

### Uploading with the API

`zenodo/upload_data.py` uploads the five files to the reserved draft, checks
each against Zenodo's MD5 (skipping files already uploaded unchanged), attaches
the draft to the community, and stops before publishing. `--metadata` also sets
the metadata from `zenodo/data.json`, keeping the reserved DOI; leave it out if
the metadata was entered in the browser. It reads a personal access token
(scope `deposit:write`) from `~/.zenodo_token`, which must be readable only by
you. On NCI Gadi, run it in the `copyq` queue, the only one with internet:

```bash
chmod 600 ~/.zenodo_token
qsub -P fy54 -q copyq -l ncpus=1,mem=4GB,walltime=02:00:00,storage=gdata/fy54 \
  -j oe -o zenodo_upload.log -- /bin/bash -c \
  'cd /g/data/fy54/analysis/tenk10k-effector-genes && python3 zenodo/upload_data.py --record 23406585 --dry-run'
```

Check the log, then run it again without `--dry-run`.

## Licences

The code is under the MIT licence (`LICENSE`). The data record is under
CC BY-SA 4.0; files from third parties keep their own terms, which
`INVENTORY.tsv` lists file by file.
