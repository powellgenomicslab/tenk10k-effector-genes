FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL org.opencontainers.image.source=https://github.com/powellgenomicslab/tenk10k-effector-genes
LABEL io.github.snakemake.conda_env_hash="c39561da9a4598840ce2d770295c38106004b392c94388ba6d26ef34f76f6622"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/figures.yaml
#   prefix: /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc
#   # Environment for the manuscript figure and table scripts in scripts/
#   # (referenced from scripts/Snakefile as "../workflow/envs/figures.yaml").
#   #
#   # R 4.4.1 matches the environment that produced the published results (see
#   # renv.yaml). Package set derived from the library() and pkg:: calls in
#   # scripts/**/*.R. Exact versions are pinned in figures.<platform>.pin.txt,
#   # which Snakemake uses automatically; regenerate them with `make lock`.
#   #
#   # Not on conda: sciscales (GitHub) and oaqc (CRAN, needed by the "backbone"
#   # graph layout) are installed by figures.post-deploy.sh.
#   # The MoMAColors and MexBrewer palettes are bundled with paletteer.
#   name: figures
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base=4.4.1
#     # core data handling
#     - r-tidyverse
#     - r-data.table
#     - r-r.utils        # fread() of .gz files
#     - r-bit64          # fread() of integer64 columns
#     - r-arrow
#     - r-readxl
#     - r-writexl
#     - r-glue
#     - r-fs
#     - r-here
#     - r-janitor
#     - r-broom
#     - r-formattable
#     - r-tidytext
#     # statistics
#     - bioconductor-qvalue
#     # plotting
#     - r-ggplot2
#     - r-patchwork
#     - r-cowplot
#     - r-ggrepel
#     - r-ggnewscale
#     - r-ggtext
#     - r-ggthemes
#     - r-ggsci
#     - r-ggprism
#     - r-ggridges
#     - r-ggforce
#     - r-ggrastr
#     - r-geomtextpath
#     - r-paletteer
#     - r-metbrewer
#     - r-rcolorbrewer
#     - r-viridislite
#     - r-scales
#     - r-labeling
#     - bioconductor-complexheatmap
#     - r-circlize
#     - r-ggraph
#     - r-tidygraph
#     - r-graphlayouts
#     - r-ragg
#     - r-svglite
#     - r-png
#     # locus plots
#     - r-locuszoomr
#     - r-ldlinkr
#     - bioconductor-ensdb.hsapiens.v86
#     - bioconductor-annotationhub
#     # Fig. 1 layout, and fitting figures to the journal's print size
#     # (scripts/util/fit_journal.R: typst, pdfinfo, pdftocairo)
#     - typst
#     - poppler
#     # compiling oaqc in figures.post-deploy.sh
#     - cxx-compiler
#     - make
RUN mkdir -p /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc
COPY workflow/envs/figures.yaml /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc/environment.yaml
COPY workflow/envs/figures.linux-64.pin.txt /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc/pin.txt
COPY workflow/envs/figures.post-deploy.sh /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc/post-deploy.sh

# Step 3: Generate conda environments

RUN conda create --yes --prefix /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc --file /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc/pin.txt && \
    conda run --prefix /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc bash /conda-envs/d0dbc9ad2051ad75ff2038ee5493aedc/post-deploy.sh && \
    conda clean --all -y
