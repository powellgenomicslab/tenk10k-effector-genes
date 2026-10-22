FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL org.opencontainers.image.source=https://github.com/powellgenomicslab/tenk10k-effector-genes
LABEL io.github.snakemake.conda_env_hash="3132e804e6f9b4f4ab62450054b74f5b7d0adc0ae4bb5ccea212a4760cbd2426"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/figures.yaml
#   prefix: /conda-envs/c4e57bd56618bc203a4f679569c6b0ed
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
#     - r-gridgraphics   # patchwork: base-graphics panels (ED Fig. 9)
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
RUN mkdir -p /conda-envs/c4e57bd56618bc203a4f679569c6b0ed
COPY workflow/envs/figures.yaml /conda-envs/c4e57bd56618bc203a4f679569c6b0ed/environment.yaml
COPY workflow/envs/figures.linux-64.pin.txt /conda-envs/c4e57bd56618bc203a4f679569c6b0ed/pin.txt
COPY workflow/envs/figures.post-deploy.sh /conda-envs/c4e57bd56618bc203a4f679569c6b0ed/post-deploy.sh

# Step 3: Generate conda environments

RUN conda create --yes --prefix /conda-envs/c4e57bd56618bc203a4f679569c6b0ed --file /conda-envs/c4e57bd56618bc203a4f679569c6b0ed/pin.txt && \
    conda run --prefix /conda-envs/c4e57bd56618bc203a4f679569c6b0ed bash /conda-envs/c4e57bd56618bc203a4f679569c6b0ed/post-deploy.sh && \
    conda clean --all -y
