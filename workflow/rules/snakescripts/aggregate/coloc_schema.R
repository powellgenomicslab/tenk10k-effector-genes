# Shared reader for per-biosample/per-pheno all_chr.coloc.tsv files.
#
# Sourced by concat_coloc_all.R and by scripts/aggregate/rebuild_coloc_parquet.R
# so both use identical normalisation (same pattern as the format_gwas
# snakescripts, which all source() a shared liftover helper).
#
# Why this is not just fread()
# ----------------------------
# run_coloc_test.R historically named the SNP count `nsnps_coloc` on its error
# path and `nsnps_coloc_tested` on its success path. bind_rows() unioned both, so
# a per-chromosome file had 12 or 13 columns with either name in position 5 or
# 13. The old concat_coloc_chr kept only the first chromosome's header, so ~80
# aggregated files carry a 12-name header over rows that are a mix of 12 and 13
# fields.
#
# On those files `fread(fill = TRUE)` commits to 12 columns from its initial row
# sample and then STOPS EARLY at the first 13-field row, silently discarding
# every remaining row - e.g. ASDC/alzheimers returns 317 of its 653 rows. That
# truncation is why the previously shipped parquet was short.
#
# So parse by line instead: split on tabs, keep the true maximum field count,
# and name any field past the header V<n>. This cannot truncate.

library(data.table)

# Canonical numeric columns; everything else stays character.
.coloc_num_cols <- c("chr", "nsnps_coloc_tested", "PP.H0.abf", "PP.H1.abf",
                     "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "top_snp_pph4")

# Number of rows whose SNP count had to be recovered from an alias column,
# accumulated across files for reporting.
.coloc_recovered <- 0L

# Parse and normalise a single file. Coalescing and type conversion happen here
# rather than after rbindlist so the combined table is never held as ~140M
# character cells.
read_one_coloc_file <- function(f) {
    lines <- readLines(f, warn = FALSE)
    if (length(lines) < 2L) return(NULL)   # header-only or empty
    hdr <- strsplit(lines[1L], "\t", fixed = TRUE)[[1L]]
    cols <- tstrsplit(lines[-1L], "\t", fixed = TRUE)
    names(cols) <- if (length(cols) <= length(hdr)) {
        hdr[seq_along(cols)]
    } else {
        c(hdr, paste0("V", seq(length(hdr) + 1L, length(cols))))
    }
    d <- as.data.table(cols)

    # Blank fields read as "" - treat as missing before coalescing.
    for (j in names(d)) d[get(j) == "", (j) := NA_character_]

    if (!"nsnps_coloc_tested" %in% names(d)) d[, nsnps_coloc_tested := NA_character_]
    aliases <- intersect(c("nsnps_coloc", "V13"), names(d))
    for (alias in aliases) {
        .coloc_recovered <<- .coloc_recovered +
            d[is.na(nsnps_coloc_tested) & !is.na(get(alias)), .N]
        d[is.na(nsnps_coloc_tested), nsnps_coloc_tested := get(alias)]
    }
    if (length(aliases)) d[, (aliases) := NULL]

    num_cols <- intersect(.coloc_num_cols, names(d))
    d[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
    if ("chr" %in% names(d)) d[, chr := as.integer(chr)]
    d
}

# Coalesce the SNP count from every place it can appear, then drop the aliases.
# Errors rather than silently emitting NA, so a new schema variant is caught.
read_coloc_all_chr <- function(files, expected_n = NULL) {
    if (!is.null(expected_n) && length(files) != expected_n) {
        stop(sprintf(
            "expected %d input files but got %d - check target_phenotypes.txt",
            expected_n, length(files)
        ))
    }
    if (!length(files)) stop("no input files given")

    .coloc_recovered <<- 0L
    combined <- rbindlist(lapply(files, read_one_coloc_file), fill = TRUE)
    if (!nrow(combined)) stop("all input files were empty")
    cat(sprintf("Recovered SNP count for %d rows from alias columns\n",
                .coloc_recovered))

    # Guard 1: no unnamed column may survive. A surviving V<n> means a schema
    # variant this function does not know about.
    stray <- grep("^V[0-9]+$", names(combined), value = TRUE)
    if (length(stray)) {
        stop("unnamed column(s) survived normalisation: ", paste(stray, collapse = ", "))
    }

    # Guard 2: every row must have a SNP count. Rows from genes that failed
    # coloc.abf are exempt - run_coloc_test.R records those with status != "ok".
    n_na <- if ("status" %in% names(combined)) {
        combined[is.na(nsnps_coloc_tested) & (is.na(status) | status == "ok"), .N]
    } else {
        combined[is.na(nsnps_coloc_tested), .N]
    }
    if (n_na > 0) {
        stop(sprintf("%d successful rows still have nsnps_coloc_tested = NA", n_na))
    }

    combined
}
