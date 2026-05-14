# ======================================================================
# STRING enrichment via REST API  (https://string-db.org/help/api/)
#
# Two modes:
#   1. SET-BASED (default) — gene files have one column (gene names).
#      Calls /api/tsv/enrichment with optional custom background.
#   2. RANK-BASED — gene files have two columns (gene <TAB> value).
#      Retrieves /api/tsv/functional_annotation for the full universe
#      and performs local Wilcoxon rank-sum (Mann-Whitney U) tests per
#      term, testing whether annotated genes rank differently from
#      un-annotated ones.
#
# Follows the same loop-over-phenotype / combined-output pattern as
# gprofiler.R and writes a TSV that is compatible with the
# aggregate_enrichment rule in enrichment.smk.
# ======================================================================
library(httr)
library(data.table)

INPUT  <- snakemake@input
OUTPUT <- snakemake@output
PARAMS <- snakemake@params
LOG    <- snakemake@log

# ------------------------------------------------------------------
# Logging boiler-plate  (identical to gprofiler.R)
# ------------------------------------------------------------------
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

log_msg("Start: STRING REST-API enrichment")
log_msg("  dir_gene_set:  ", INPUT$dir_gene_set)
log_msg("  gene_universe: ", INPUT$gene_universe)
log_msg("  output:        ", OUTPUT$enrich)

# ------------------------------------------------------------------
# Parameters (with sensible defaults)
# ------------------------------------------------------------------
species        <- if (!is.null(PARAMS$species))        PARAMS$species        else 9606
string_version <- if (!is.null(PARAMS$string_version)) PARAMS$string_version else "12.0"
min_genes      <- if (!is.null(PARAMS$min_genes))      PARAMS$min_genes      else 5
fdr_threshold  <- if (!is.null(PARAMS$fdr_threshold))  PARAMS$fdr_threshold  else 0.05
caller_id      <- if (!is.null(PARAMS$caller_identity)) PARAMS$caller_identity else "tenk10k_smr"

# Build the base URL pinned to the requested STRING version
string_api_url <- paste0("https://version-", gsub("\\.", "-", string_version),
                         ".string-db.org/api")

log_msg("  species:        ", species)
log_msg("  string_version: ", string_version)
log_msg("  string_api_url: ", string_api_url)
log_msg("  min_genes:      ", min_genes)
log_msg("  fdr_threshold:  ", fdr_threshold)
log_msg("  caller_identity:", caller_id)

# ------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------
t0 <- Sys.time()

files <- list.files(INPUT$dir_gene_set, pattern = "\\.txt$", full.names = TRUE)
names(files) <- tools::file_path_sans_ext(basename(files))

gene_universe <- readLines(INPUT$gene_universe)

log_msg("Found ", length(files), " phenotype file(s) in: ", INPUT$dir_gene_set)
log_msg("Gene universe size: ", length(gene_universe))

# ===================================================================
# Helper: POST to STRING API with polite 1-second sleep
# ===================================================================
string_post <- function(method, params, format = "tsv") {
    url <- paste(string_api_url, format, method, sep = "/")
    Sys.sleep(1)                       # STRING asks for >= 1 s between calls
    resp <- httr::POST(url, body = params, encode = "form")
    httr::stop_for_status(resp)
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (nchar(trimws(txt)) == 0) return(data.table::data.table())
    data.table::fread(text = txt, sep = "\t", header = TRUE)
}

# ===================================================================
# Helper: map gene names → STRING IDs   (batched, max ~2000 / call)
# ===================================================================
map_identifiers <- function(genes, batch_size = 2000) {
    all_mapped <- list()
    for (start in seq(1, length(genes), by = batch_size)) {
        batch <- genes[start:min(start + batch_size - 1, length(genes))]
        res <- tryCatch(
            string_post("get_string_ids", list(
                identifiers   = paste(batch, collapse = "\r"),
                species       = species,
                limit         = 1,
                echo_query    = 1,
                caller_identity = caller_id
            )),
            error = function(e) {
                log_msg("  WARN map_identifiers batch ", start, ": ", conditionMessage(e))
                data.table::data.table()
            }
        )
        if (nrow(res) > 0) all_mapped <- c(all_mapped, list(res))
    }
    data.table::rbindlist(all_mapped, fill = TRUE)
}

# ===================================================================
# Map the background / universe  (once)
# ===================================================================
log_msg("Mapping gene universe to STRING IDs...")
universe_mapped <- map_identifiers(gene_universe)
if (nrow(universe_mapped) == 0) stop("No universe genes could be mapped to STRING IDs.")
background_ids <- unique(universe_mapped$stringId)
log_msg("  Mapped: ", length(background_ids), " / ", length(gene_universe), " universe genes")

# ===================================================================
# Read a gene-set file.  Returns list(genes, values, ranked)
# If file has 2 tab/space columns → ranked mode
# ===================================================================
read_gene_file <- function(path) {
    lines <- readLines(path)
    lines <- lines[nzchar(trimws(lines))]
    # check if two-column (gene <TAB/space> value)
    first <- strsplit(lines[1], "\\s+")[[1]]
    if (length(first) >= 2 && !is.na(suppressWarnings(as.numeric(first[2])))) {
        parsed <- data.table::fread(text = paste(lines, collapse = "\n"),
                                    header = FALSE, col.names = c("gene", "value"))
        parsed <- parsed[order(-abs(value))]      # sort by |value| descending
        list(genes = parsed$gene, values = parsed$value, ranked = TRUE)
    } else {
        list(genes = lines, values = NULL, ranked = FALSE)
    }
}

# ===================================================================
# MODE 1: Set-based enrichment via /api/tsv/enrichment
# ===================================================================
run_set_enrichment <- function(query_ids) {
    params <- list(
        identifiers                  = paste(query_ids, collapse = "\r"),
        species                      = species,
        background_string_identifiers = paste(background_ids, collapse = "\r"),
        caller_identity              = caller_id
    )
    res <- tryCatch(
        string_post("enrichment", params),
        error = function(e) {
            log_msg("  ERROR enrichment API: ", conditionMessage(e))
            data.table::data.table()
        }
    )
    if (nrow(res) == 0) return(NULL)

    # Rename columns → gprofiler-compatible schema
    col_map <- c(
        source             = "category",
        term_id            = "term",
        term_name          = "description",
        p_value            = "p_value",
        fdr                = "fdr",
        term_size          = "number_of_genes_in_background",
        intersection_size  = "number_of_genes",
        intersection_genes = "inputGenes",
        query_gene_names   = "preferredNames"
    )
    for (new_nm in names(col_map)) {
        old_nm <- col_map[[new_nm]]
        if (old_nm %in% names(res)) data.table::setnames(res, old_nm, new_nm)
    }
    res[, query_size  := length(query_ids)]
    res[, significant := fdr < fdr_threshold]
    res[, enrichment_method := "set"]
    res <- res[significant == TRUE]
    if (nrow(res) == 0) return(NULL)
    res
}

# ===================================================================
# PPI enrichment via /api/tsv/ppi_enrichment
# Tests whether the query network has more interactions than expected
# ===================================================================
run_ppi_enrichment <- function(query_ids) {
    params <- list(
        identifiers                   = paste(query_ids, collapse = "\r"),
        species                       = species,
        background_string_identifiers = paste(background_ids, collapse = "\r"),
        caller_identity               = caller_id
    )
    res <- tryCatch(
        string_post("ppi_enrichment", params),
        error = function(e) {
            log_msg("  ERROR ppi_enrichment API: ", conditionMessage(e))
            data.table::data.table()
        }
    )
    if (nrow(res) == 0) return(NULL)

    # The endpoint returns a single row with:
    #   number_of_nodes, number_of_edges, average_node_degree,
    #   local_clustering_coefficient, expected_number_of_edges, p_value
    dt <- data.table::data.table(
        source                       = "PPI",
        term_id                      = "PPI_enrichment",
        term_name                    = "Protein-protein interaction enrichment",
        p_value                      = res$p_value,
        fdr                          = res$p_value,   # single test, no correction needed
        term_size                    = NA_integer_,
        intersection_size            = res$number_of_edges,
        intersection_genes           = NA_character_,
        query_gene_names             = NA_character_,
        query_size                   = res$number_of_nodes,
        significant                  = res$p_value < fdr_threshold,
        enrichment_method            = "ppi",
        ppi_expected_edges           = res$expected_number_of_edges,
        ppi_avg_node_degree          = res$average_node_degree,
        ppi_local_clustering_coeff   = res$local_clustering_coefficient
    )
    dt
}

# ===================================================================
# MODE 2: Rank-based enrichment via functional_annotation + local
#          Wilcoxon rank-sum / Mann-Whitney U test
# ===================================================================

# Cache: retrieve functional annotations for the universe (once)
anno_cache <- NULL

get_annotations <- function() {
    if (!is.null(anno_cache)) return(anno_cache)
    log_msg("Retrieving functional annotations for universe from STRING...")
    params <- list(
        identifiers     = paste(background_ids, collapse = "\r"),
        species         = species,
        caller_identity = caller_id
    )
    # functional_annotation may be large; fetch in batches if needed
    anno <- tryCatch(
        string_post("functional_annotation", params),
        error = function(e) {
            log_msg("  ERROR functional_annotation API: ", conditionMessage(e))
            data.table::data.table()
        }
    )
    anno_cache <<- anno
    log_msg("  Retrieved ", nrow(anno), " annotation rows for ", 
            length(unique(anno$inputGenes)), " unique gene entries")
    anno
}

run_rank_enrichment <- function(genes, values) {
    anno <- get_annotations()
    if (nrow(anno) == 0) {
        log_msg("  No annotations retrieved — cannot do rank-based enrichment")
        return(NULL)
    }

    # Build a gene → value lookup (using preferredNames from mapping)
    query_mapped <- map_identifiers(genes)
    if (nrow(query_mapped) == 0) return(NULL)

    # Merge values onto mapped genes
    val_dt <- data.table::data.table(gene = genes, value = values)
    query_mapped[, queryItem := as.character(queryItem)]
    val_dt[, gene := as.character(gene)]
    merged <- merge(query_mapped, val_dt, by.x = "queryItem", by.y = "gene", all.x = TRUE)
    merged <- merged[!is.na(value)]
    if (nrow(merged) == 0) return(NULL)

    # Create lookup: preferredName → value
    gene_vals <- setNames(merged$value, merged$preferredName)

    # Explode annotations:  each row's inputGenes / preferredNames is comma-separated
    # We need to test per (category, term): are values of annotated genes different?
    anno_long <- anno[, {
        prefs <- strsplit(as.character(preferredNames), ",")[[1]]
        list(preferredName = trimws(prefs))
    }, by = .(category, term, description)]

    # For each term, do a Wilcoxon rank-sum test
    terms <- unique(anno_long[, .(category, term, description)])
    results <- list()

    for (j in seq_len(nrow(terms))) {
        cat_j   <- terms$category[j]
        term_j  <- terms$term[j]
        desc_j  <- terms$description[j]

        annotated_genes <- anno_long[category == cat_j & term == term_j, preferredName]
        in_set  <- intersect(annotated_genes, names(gene_vals))
        out_set <- setdiff(names(gene_vals), annotated_genes)

        if (length(in_set) < 2 || length(out_set) < 2) next

        wt <- suppressWarnings(wilcox.test(
            gene_vals[in_set], gene_vals[out_set],
            alternative = "two.sided", exact = FALSE
        ))
        results[[length(results) + 1]] <- data.table::data.table(
            source             = cat_j,
            term_id            = term_j,
            term_name          = desc_j,
            p_value            = wt$p.value,
            term_size          = length(annotated_genes),
            intersection_size  = length(in_set),
            intersection_genes = paste(in_set, collapse = ","),
            query_gene_names   = paste(in_set, collapse = ","),
            query_size         = length(gene_vals),
            statistic          = as.numeric(wt$statistic)
        )
    }

    if (length(results) == 0) return(NULL)
    dt <- data.table::rbindlist(results)

    # FDR correction
    dt[, fdr := p.adjust(p_value, method = "BH")]
    dt[, significant := fdr < fdr_threshold]
    dt[, enrichment_method := "rank"]
    dt <- dt[significant == TRUE]
    dt <- dt[order(fdr)]
    if (nrow(dt) == 0) return(NULL)
    dt
}

# ===================================================================
# Loop over phenotypes
# ===================================================================
all_results <- list()

for (i in seq_along(files)) {
    pheno   <- names(files)[i]
    file    <- files[i]
    t_pheno <- Sys.time()

    log_msg("[", i, "/", length(files), "] ", pheno,
            ": loading gene set from ", file)

    gs <- read_gene_file(file)
    log_msg("[", i, "/", length(files), "] ", pheno,
            ": ", length(gs$genes), " genes, ranked=", gs$ranked)

    if (length(gs$genes) < min_genes) {
        log_msg("[", i, "/", length(files), "] ", pheno,
                ": skipping (<", min_genes, " genes)")
        next
    }

    # -- Map query genes to STRING IDs (needed for both modes + PPI) --
    log_msg("[", i, "/", length(files), "] ", pheno,
            ": mapping query genes to STRING IDs...")
    query_mapped <- map_identifiers(gs$genes)
    query_ids    <- unique(query_mapped$stringId)
    log_msg("[", i, "/", length(files), "] ", pheno,
            ": mapped ", length(query_ids), "/", length(gs$genes), " genes")

    if (length(query_ids) == 0) {
        log_msg("[", i, "/", length(files), "] ", pheno, ": no genes mapped")
        next
    }

    if (gs$ranked) {
        # ----- Rank-based enrichment -----
        log_msg("[", i, "/", length(files), "] ", pheno,
                ": running rank-based enrichment (Wilcoxon)...")
        dt <- tryCatch(
            run_rank_enrichment(gs$genes, gs$values),
            error = function(e) {
                log_msg("ERROR [", pheno, " rank]: ", conditionMessage(e))
                NULL
            }
        )
    } else {
        # ----- Set-based enrichment via STRING API -----
        log_msg("[", i, "/", length(files), "] ", pheno,
                ": calling STRING enrichment API...")
        dt <- tryCatch(
            run_set_enrichment(query_ids),
            error = function(e) {
                log_msg("ERROR [", pheno, " set]: ", conditionMessage(e))
                NULL
            }
        )
    }

    # ----- PPI enrichment (always, regardless of mode) -----
    log_msg("[", i, "/", length(files), "] ", pheno,
            ": calling STRING PPI enrichment API...")
    dt_ppi <- tryCatch(
        run_ppi_enrichment(query_ids),
        error = function(e) {
            log_msg("ERROR [", pheno, " ppi]: ", conditionMessage(e))
            NULL
        }
    )
    if (!is.null(dt_ppi) && nrow(dt_ppi) > 0) {
        log_msg("[", i, "/", length(files), "] ", pheno,
                ": PPI p=", dt_ppi$p_value,
                ", edges=", dt_ppi$intersection_size,
                "/expected=", dt_ppi$ppi_expected_edges)
    }

    # Combine functional + PPI results for this phenotype
    dt <- data.table::rbindlist(list(dt, dt_ppi), fill = TRUE)

    if (is.null(dt) || nrow(dt) == 0) {
        log_msg("[", i, "/", length(files), "] ", pheno,
                ": no significant enrichment results")
    } else {
        n_terms <- nrow(dt)
        n_sig   <- sum(dt$significant, na.rm = TRUE)
        dt[, phenotype := pheno]
        all_results <- c(all_results, list(dt))
        elapsed_pheno <- round(difftime(Sys.time(), t_pheno, units = "secs"), 1)
        log_msg("[", i, "/", length(files), "] ", pheno,
                ": ", n_terms, " terms, ", n_sig,
                " significant (", elapsed_pheno, "s)")
    }
}

# ------------------------------------------------------------------
# Write combined output
# ------------------------------------------------------------------
combined <- data.table::rbindlist(all_results, fill = TRUE)
n_phenos_with_results <- length(all_results)

if (nrow(combined) == 0) {
    log_msg("WARNING: no enrichment results across all phenotypes — writing empty file")
    data.table::fwrite(data.table::data.table(), OUTPUT$enrich, sep = "\t")
} else {
    # Flatten any remaining list columns (safety net)
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
