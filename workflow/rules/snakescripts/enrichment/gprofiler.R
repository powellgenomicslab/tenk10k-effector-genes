# enrichment analysis using gProfiler gene sets
library(gprofiler2)
library(data.table)

INPUT  <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params
LOG    <- snakemake@log

# Redirect all stdout and stderr to the Snakemake log file
log_file <- LOG[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con  <- file(log_file, open = "wt")
sink(log_con,              type = "output")
sink(log_con, append = TRUE, type = "message")
on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
}, add = TRUE)

log_msg <- function(...) {
    message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

log_msg("Start: gProfiler enrichment")
log_msg("  directory of gene_set:     ", INPUT$dir_gene_set)
log_msg("  universe:     ", INPUT$gene_universe)
log_msg("  sources:      ", paste(PARAMS$sources, collapse = ", "))
log_msg("  output:       ", OUTPUT$enrich)

# ------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------
t0 <- Sys.time()

files <- list.files(INPUT$dir_gene_set, pattern = "\\.txt$", full.names = TRUE)
names(files) <- tools::file_path_sans_ext(basename(files))

gene_universe <- readLines(INPUT$gene_universe)

# ------------------------------------------------------------------
# Run gProfiler
# ------------------------------------------------------------------
enrich_gprofiler <- function(gene_set,
               custom_bg    = gene_universe,
               sources      = PARAMS$sources,
               organism     = "hsapiens",
               ordered_query = as.logical(PARAMS$ordered_query),
               multi_query  = as.logical(PARAMS$multi_query),
               highlight    = as.logical(PARAMS$highlight),
               domain_scope = PARAMS$domain_scope,
               significant = PARAMS$significant,
               ...) {
    gost(
        query         = gene_set,
        organism      = organism,
        ordered_query = ordered_query,
        multi_query   = multi_query,
        sources       = sources,
        highlight     = highlight,
        custom_bg     = custom_bg,
        domain_scope  = domain_scope,
        significant   = significant,
        ...
    )
}

log_msg("Found ", length(files), " phenotype file(s) in: ", INPUT$dir_gene_set)

# ------------------------------------------------------------------
# Loop over phenotypes
# ------------------------------------------------------------------
all_results <- list()

for (i in seq_along(files)) {
    pheno   <- names(files)[i]
    file    <- files[i]
    t_pheno <- Sys.time()

    log_msg("[", i, "/", length(files), "] ", pheno, ": loading gene set from ", file)
    gene_set <- readLines(file)
    log_msg("[", i, "/", length(files), "] ", pheno, ": ", length(gene_set), " genes")

    log_msg("[", i, "/", length(files), "] ", pheno, ": running gProfiler...")
    res <- tryCatch(
        enrich_gprofiler(gene_set),
        error = function(e) {
            log_msg("ERROR [", pheno, "]: ", conditionMessage(e))
            NULL
        }
    )

    if (!is.null(PARAMS$custom_gmt) && length(PARAMS$custom_gmt) > 0) {
        for (gmt_key in names(PARAMS$custom_gmt)) {
            gmt_value <- PARAMS$custom_gmt[[gmt_key]]
            log_msg("[", i, "/", length(files), "] ", pheno, ": running gProfiler with custom GMT: ", gmt_key)
            res_custom <- tryCatch(
                enrich_gprofiler(gene_set, organism = gmt_value),
                error = function(e) {
                    log_msg("ERROR [", pheno, ", GMT: ", gmt_key, "]: ", conditionMessage(e))
                    NULL
                }
            )
            
            if (!is.null(res_custom) && !is.null(res_custom$result) && nrow(res_custom$result) > 0) {
                n_terms <- nrow(res_custom$result)
                n_sig   <- sum(res_custom$result$significant, na.rm = TRUE)
                dt      <- data.table::as.data.table(res_custom$result)
                dt[, phenotype := pheno]
                dt[, custom_gmt := gmt_key]
                all_results <- c(all_results, list(dt))
                log_msg("[", i, "/", length(files), "] ", pheno, " (", gmt_key, "): ",
                        n_terms, " terms, ", n_sig, " significant")
            }
        }
    }

    if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
        log_msg("[", i, "/", length(files), "] ", pheno, ": no enrichment results")
    } else {
        n_terms <- nrow(res$result)
        n_sig   <- sum(res$result$significant, na.rm = TRUE)
        dt      <- data.table::as.data.table(res$result)
        dt[, phenotype := pheno]
        all_results <- c(all_results, list(dt))
        elapsed_pheno <- round(difftime(Sys.time(), t_pheno, units = "secs"), 1)
        log_msg("[", i, "/", length(files), "] ", pheno, ": ",
                n_terms, " terms, ", n_sig, " significant (", elapsed_pheno, "s)")
    }
}

# ------------------------------------------------------------------
# Write combined output
# ------------------------------------------------------------------
combined <- data.table::rbindlist(all_results, fill = TRUE)
n_phenos_with_results <- sum(!sapply(all_results, is.null))

if (nrow(combined) == 0) {
    log_msg("WARNING: no enrichment results across all phenotypes — writing empty output")
    data.table::fwrite(data.table::data.table(), OUTPUT$enrich, sep = "\t")
} else {
    # Unnest list columns
    list_cols <- names(combined)[sapply(combined, is.list)]
    for (col in list_cols) {
        combined[[col]] <- sapply(combined[[col]], function(x) {
            if (length(x) == 0) NA_character_ else paste(as.character(x), collapse = ";")
        })
    }
    
    log_msg("Writing combined results: ", nrow(combined), " rows from ",
            n_phenos_with_results, "/", length(files), " phenotypes")
    data.table::fwrite(combined, OUTPUT$enrich, sep = "\t")
}

elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
log_msg("Done: ", OUTPUT$enrich, " (", elapsed, "s)")
