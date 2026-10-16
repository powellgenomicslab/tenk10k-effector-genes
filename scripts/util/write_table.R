# misc helper function

library(dplyr)

write_table <- function(df, source_tbl_name,
                        tbl_order = NULL,
                        target_tbl_name = source_tbl_name,
                        target_dir = "tables",
                        cols = NULL,
                        label_xlsx = "metadata/table_column_names.xlsx") {
    if (!is.null(cols)) {
        df <- df %>% select(all_of(cols))
    }
    # manuscript Supplementary Table N is written to tables/tab-sN.tsv
    if (!is.null(tbl_order)) {
        target_tbl_name <- paste0("tab-s", tbl_order)
    }
    labels <- readxl::read_excel(label_xlsx, sheet = source_tbl_name) %>%
        select(label, name) %>%
        deframe()
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    df %>%
        select(!!!labels) %>%
        data.table::fwrite(file.path(target_dir, paste0(target_tbl_name, ".tsv")), sep = "\t", na = "NA", quote = FALSE)
}
