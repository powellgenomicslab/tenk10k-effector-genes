library(schard)
library(tidyverse)
library(data.table)
library(scPred)
# library(BPCells)
library(SeuratObject)
library(Seurat)

# interactive
STUDY <- snakemake@wildcards[["study"]]
ANNOT <- snakemake@wildcards[["annot"]]
H5AD <- snakemake@input[["h5ad"]]
# DIR_COUNTS <- snakemake@input[["dir_counts"]]
REF_SCPRED <- snakemake@input[["ref_scpred"]]
PARAMS <- snakemake@params
OUTPUT <- snakemake@output

# STUDY <- "immune_atlas"
# ANNOT <- "scpred"
# H5AD <- glue::glue("resources/scdrs/h5ad/{STUDY}.prep.h5ad")
# DIR_COUNTS <- glue::glue("resources/scdrs/counts/{STUDY}/")
# REF_SCPRED <- glue::glue("resources/scpred/{STUDY}.ref.rds")
# PARAMS <- list(future_globals_maxSize = 640000 * 1024^2,
#                scpred_max_iter_harmony = 20,
#                scpred_threshold = 0.5)

# Load data use BPCells functions
# fs::dir_create(DIR_COUNTS)
# query <- open_matrix_annquery_hdf5(path = H5AD)



# # Write the matrix to a directory
# write_matrix_dir(
#   mat = query,
#   dir = DIR_COUNTS,
#   overwrite = TRUE
# )

# # Now that we have the matrix on disk, we can load it
# mat <- open_matrix_dir(dir = DIR_COUNTS)
# data <- CreateSeuratObject(counts = mat)

# Load query with schard (if small enough)
query <- schard::h5ad2seurat(H5AD)

# map to HGNC symbols
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
gene_map <- mapIds(EnsDb.Hsapiens.v86, 
                   keys = rownames(query), 
                   keytype = "GENEID", 
                   column = "SYMBOL", 
                   multiVals = "first")

# Replace genes that do not have a mapping with original rownames
gene_map[is.na(gene_map)] <- rownames(query)[is.na(gene_map)]

# ensure unique gene names 
gene_map <- make.unique(gene_map)

# Update rownames of Seurat object
mat_counts <- LayerData(query, assay = "RNA", layer = "counts")
mat_data <- LayerData(query, assay = "RNA", layer = "data")
rownames(mat_counts) <- gene_map
rownames(mat_data) <- gene_map

# Create Assay with counts and data
query <- CreateSeuratObject(
  counts = mat_counts,
  data = mat_data
)

# Query data preprocessing
query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^MT-")
query <- NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000)
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(query)

query <- ScaleData(query)

# Unable to run these two lines due to memory issues
query <- SCTransform(query)
DefaultAssay(query) <- "SCT"

query <- RunPCA(query, npcs = 10) #, assay = "RNA"
query <- FindNeighbors(query,dims = 1:10, k.param = 10)
query <- FindClusters(query)

# WARNING
# Warning message:
# UNRELIABLE VALUE: One of the _future.apply_ iterations (_future_lapply-1_) unexpectedly generated random numbers without declaring so. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'future.seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via a parallel RNG method. To disable this check, use 'future.seed = NULL', or set option 'future.rng.onMisuse' to "ignore". 
query <- RunUMAP(query,dims = 1:10)

# save data
# saveRDS(query, glue::glue("results/scdrs/annot_sample/{STUDY}.query.rds"))


# plan("multicore", workers = as.numeric(Sys.getenv("NCPUS")))
# options(future.globals.maxSize = PARAMS$future_globals_maxSize)

# # prepare reference for scPred
# ref <- readRDS(REF_SCPRED)

# # hierarchical scPred annotation
# library(data.tree)
# library(HierscPred)
# query <- HierscPred::predictTree(
#     ref,
#     newData = query,
#     threshold = PARAMS$scpred_threshold,
#     max.iter.harmony = PARAMS$scpred_max_iter_harmony
# )

query <- scPredict(query, ref, max.iter.harmony = PARAMS$scpred_max_iter_harmony, threshold = PARAMS$scpred_threshold)

# OUTPUT <- list(
#   raw = glue::glue("results/scdrs/annot_sample/{STUDY}.{ANNOT}.raw.csv"),
#   format = glue::glue("results/scdrs/annot_sample/{STUDY}.{ANNOT}.format.csv")
# )


fs::dir_create(dirname(OUTPUT$raw))

df <- query@meta.data %>%
  rownames_to_column(var = "cell_id")

fwrite(df, OUTPUT$raw, row.names = FALSE)

df_format <- df %>%
  select(cell_id, cell_type = scpred_prediction) |> 
  mutate(cell_type = str_replace_all(cell_type, " ", "_"))

fwrite(df_format, OUTPUT$format, row.names = FALSE)