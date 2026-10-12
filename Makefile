# Reproduce the manuscript figures and tables.
#
#   make fig-m3        Main Figure 3              -> figures/fig-m3.{pdf,png}
#   make fig-e3        Extended Data Figure 3     -> figures/fig-e3.{pdf,png}
#   make fig-s3        Supplementary Figure 3     -> figures/fig-s3.{pdf,png}
#   make tab-s3        Supplementary Table 3      -> tables/tab-s3.tsv
#   make figs-m | figs-e | figs-s | tabs-s | all
#   make list          every target with its caption
#
# Inputs: `make data` downloads the public summary-statistics bundle, or
# `make link-data DATA_DIR=/path/to/tenk10k-causal` links the files listed in
# config/figure_inputs.txt from an existing results tree.
#
# Software: SDM=conda (default) builds the pinned conda environment
# (workflow/envs/figures.<platform>.pin.txt); once the container image is
# published, SDM="conda apptainer" runs inside it; SDM= uses software on PATH.

SNAKEMAKE ?= snakemake
CORES     ?= 1
SDM       ?= conda
ARGS      ?=
SMK = $(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none \
      --cores $(CORES) $(if $(SDM),--sdm $(SDM)) $(ARGS) --

.PHONY: help list all figs-m figs-e figs-s tabs-s unnumbered data link-data \
        manifest env lock dockerfile dry-run

help:
	@sed -n '2,17p' Makefile | sed 's/^# \{0,1\}//'

list:
	@$(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none --cores 1 \
	  --quiet all -- targets

all figs-m figs-e figs-s tabs-s unnumbered:
	$(SMK) $(subst -,_,$@)

fig-% tab-%:
	$(SMK) $(subst -,_,$@)

dry-run:
	$(MAKE) all ARGS=--dry-run

data:
	$(SMK) download_inputs

DATA_DIR ?=
link-data: config/figure_inputs.txt
	@test -n "$(DATA_DIR)" || { echo "usage: make link-data DATA_DIR=/path/to/data/root"; exit 1; }
	bash scripts/0-preprocess/link_figure_inputs.sh "$(DATA_DIR)"

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

dockerfile:
	$(SNAKEMAKE) --snakefile scripts/Snakefile --workflow-profile none \
	  --containerize dockerfile -- containerize > Dockerfile
	python3 workflow/envs/pin_dockerfile.py Dockerfile
