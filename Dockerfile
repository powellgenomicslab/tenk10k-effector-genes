FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL org.opencontainers.image.source=https://github.com/powellgenomicslab/tenk10k-effector-genes
LABEL io.github.snakemake.conda_env_hash="14336bddede0507f63ea45e8de10b4e392b7100488ef5b0e0c3dfb89cbbe1a09"
# Conda environment:
#   source: workflow/envs/figures.yaml
#   prefix: /conda-envs/0e5d19ca8ec224084e72c6f4857c698c
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
RUN mkdir -p /conda-envs/0e5d19ca8ec224084e72c6f4857c698c
COPY workflow/envs/figures.yaml /conda-envs/0e5d19ca8ec224084e72c6f4857c698c/environment.yaml
COPY workflow/envs/figures.linux-64.pin.txt /conda-envs/0e5d19ca8ec224084e72c6f4857c698c/pin.txt
COPY workflow/envs/figures.post-deploy.sh /conda-envs/0e5d19ca8ec224084e72c6f4857c698c/post-deploy.sh
# Conda environment:
#   source: workflow/envs/scverse.yaml
#   prefix: /conda-envs/4c190834c12b1a68e32971e995c9caf0
#   # Single-cell Python environment for the scDRS rules (`conda: "scverse"`).
#   name: scverse
#   channels:
#     - conda-forge
#   dependencies:
#     - python=3.11
#     - scanpy
#     - anndata
#     - python-igraph
#     - leidenalg
#     - pyarrow
#     - numpy
#     - pandas
#     - pyyaml
#     - git
#     - pip
#     - pip:
#       # scDRS is not packaged for conda, and 1.0.3 is not on PyPI (which stops at
#       # 1.0.2): install the v1.0.3 release tag from GitHub, pinned to its commit
#       - scdrs @ git+https://github.com/martinjzhang/scDRS@bc1e77b58ec53296857701018a86d4494e23806b
RUN mkdir -p /conda-envs/4c190834c12b1a68e32971e995c9caf0
COPY workflow/envs/scverse.yaml /conda-envs/4c190834c12b1a68e32971e995c9caf0/environment.yaml

RUN conda create --yes --prefix /conda-envs/0e5d19ca8ec224084e72c6f4857c698c --file /conda-envs/0e5d19ca8ec224084e72c6f4857c698c/pin.txt && \
    conda run --prefix /conda-envs/0e5d19ca8ec224084e72c6f4857c698c bash /conda-envs/0e5d19ca8ec224084e72c6f4857c698c/post-deploy.sh && \
    conda env create --prefix /conda-envs/4c190834c12b1a68e32971e995c9caf0 --file /conda-envs/4c190834c12b1a68e32971e995c9caf0/environment.yaml && \
    conda clean --all -y
