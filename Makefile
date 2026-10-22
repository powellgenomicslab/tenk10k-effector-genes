# Reproduce the manuscript figures and tables.
#
#   make fig-m3        Main Figure 3              -> figures/fig-m3.png
#   make fig-e3        Extended Data Figure 3     -> figures/fig-e3.png
#   make fig-s3        Supplementary Figure 3     -> figures/fig-s3.png
#   make tab-s3        Supplementary Table 3      -> tables/tab-s3.tsv
#   make figs-m | figs-e | figs-s | tabs-s | all
#   make figures/panels/<folder>/<file>   one output file, e.g. a panel
#   make list          every target with its caption
#   make pull          download the container image only (needs internet)
#   make inventory     list the bundle files with checksums -> INVENTORY.tsv
#   make bundle        pack the bundle for Zenodo (BUNDLE=out.tar.gz)
#
# Figures are PNG; FORMATS=pdf or FORMATS=both for PDF or both. CORES=N runs N
# scripts at once (about 16 GB of memory each).
#
# Inputs: `make data` downloads the public summary-statistics bundle, or
# `make link-data DATA_DIR=/path/to/tenk10k-causal` links the files listed in
# config/figure_inputs.txt from an existing results tree.
#
# Software: SDM=conda (default) builds the pinned conda environment
# (workflow/envs/figures.<platform>.pin.txt); SDM="conda apptainer" runs
# inside the container image; SDM= uses software on PATH.

# Run from the physical path of the checkout. Apptainer starts each step in $PWD,
# and Snakemake replaces the home directory inside the container with the
# checkout; a checkout reached through a symlink in the home directory would
# then not exist inside the container.
export PWD := $(CURDIR)

SNAKEMAKE ?= snakemake
CORES     ?= 1
SDM       ?= conda
ARGS      ?=
FORMATS   ?=
SMK = $(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none \
      --cores $(CORES) $(if $(SDM),--sdm $(SDM)) $(if $(FORMATS),--config formats=$(FORMATS)) $(ARGS) --

.PHONY: help list all figs-m figs-e figs-s tabs-s data link-data \
        manifest inventory bundle env lock dockerfile dry-run pull FORCE

help:
	@sed -n '2,24p' Makefile | sed 's/^# \{0,1\}//'

list:
	@$(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none --cores 1 \
	  --quiet all -- targets

all figs-m figs-e figs-s tabs-s:
	$(SMK) $(subst -,_,$@)

fig-% tab-%:
	$(SMK) $(subst -,_,$@)

# Any output file by path, e.g. a panel: make figures/panels/1-overview/panel_b.svg
# FORCE makes make pass the request on even if the file exists; Snakemake then
# decides whether it is up to date.
figures/% tables/%: FORCE
	$(SMK) $@

FORCE:

dry-run:
	$(MAKE) all ARGS=--dry-run

# Download the container image to the file Snakemake looks for (the MD5 of its
# URL, under APPTAINER_CACHEDIR or .snakemake/singularity), so that builds can
# then run on nodes without internet access.
IMAGE = $(shell sed -n 's/^container: *"\(.*\)"/\1/p' config/figures.yaml)
pull:
	@dir="$${APPTAINER_CACHEDIR:-.snakemake/singularity}"; \
	name="$$(printf %s '$(IMAGE)' | md5sum | cut -c1-32).simg"; \
	mkdir -p "$$dir"; \
	if [ -e "$$dir/$$name" ]; then echo "$$dir/$$name is already present"; \
	else cd "$$dir" && apptainer pull --name "$$name" '$(IMAGE)'; fi

# Download and verify the data bundle. Runs the script directly, not through
# Snakemake, so that it never runs inside the container.
BUNDLE_URL = $(shell sed -n 's/^bundle_url: *"\(.*\)"/\1/p' config/figures.yaml)
data:
	bash scripts/0-preprocess/download_figure_bundle.sh '$(BUNDLE_URL)'

DATA_DIR ?=
link-data: config/figure_inputs.txt
	@test -n "$(DATA_DIR)" || { echo "usage: make link-data DATA_DIR=/path/to/data/root"; exit 1; }
	bash scripts/0-preprocess/link_figure_inputs.sh "$(DATA_DIR)"

inventory:
	$(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none --cores 1 \
	  --quiet all -- inventory

BUNDLE ?= tenk10k_figure_inputs.tar.gz
BUNDLE_ARGS ?=
bundle: inventory
	bash scripts/0-preprocess/make_figure_bundle.sh config/figure_inputs.txt $(BUNDLE) $(BUNDLE_ARGS)

manifest config/figure_inputs.txt:
	$(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none --cores 1 \
	  --quiet all -- manifest

env:
	conda create --yes --prefix .snakemake/figures-env \
	  --file workflow/envs/figures.$(shell conda info --json | python3 -c \
	  "import json,sys; print(json.load(sys.stdin)['platform'])").pin.txt
	conda run --prefix .snakemake/figures-env bash workflow/envs/figures.post-deploy.sh

lock:
	conda-lock lock --mamba -f workflow/envs/figures.yaml -p linux-64 -p osx-arm64 \
	  --kind explicit --filename-template 'workflow/envs/figures.{platform}.pin.txt'

# written to Dockerfile.new first, so that a failure leaves Dockerfile intact
dockerfile:
	$(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none \
	  --containerize -- containerize > Dockerfile.new
	python3 workflow/envs/pin_dockerfile.py Dockerfile.new
	mv Dockerfile.new Dockerfile
