args <- commandArgs(trailingOnly = TRUE)
#args <- c("breast_cimba", "pDC")
trait <- args[1]
cell_type <- args[2]
input_dir <- args[3]
output_filepath <- args[4]

library(data.table)
setDTthreads(7)

# Get all files from MR-link-2 outputs for the trait and cell type
input_filepaths <- Sys.glob(file.path(input_dir, "chr*", "*.txt"))

# Process inputs - extract names from filenames
input_dt_list <- lapply(input_filepaths, function(f) if(file.exists(f)) fread(f) else NULL)
input_dt <- rbindlist(input_dt_list)

# Get information from exposure and outcome file information
input_dt[, c("probeID", "phenotype", "biosample") := list(sub(".*(ENSG[0-9]+).*", "\\1", exposure_file),sub("\\.[^.]+$", "", basename(outcome_file)), cell_type)]
input_dt[, c("exposure_file", "outcome_file") := list(NULL, NULL)]
input_dt[, chromosome := tstrsplit(region, ":", fixed = TRUE, keep = 1L)]

input_dt <- input_dt[, .(biosample, phenotype, probeID, chromosome, region, 
    var_explained, m_snps_overlap, alpha, `se(alpha)`, `p(alpha)`, sigma_y, `se(sigma_y)`, 
    `p(sigma_y)`, sigma_x, function_time)]
input_dt[, chromosome := as.integer(chromosome)]
setorder(input_dt, chromosome)
fwrite(input_dt, output_filepath)