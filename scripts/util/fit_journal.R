# Fit a saved PDF figure to the journal's maximum print size.
#
# fit_journal_pdf() scales the whole page down uniformly (text, lines and
# points together, aspect ratio kept) so it fits within max_w_mm x max_h_mm.
# The PDF is re-wrapped with Typst, which embeds it as vectors and keeps text
# editable. Pages that already fit are left unchanged. If `png` is given, the
# PNG is re-rendered from the fitted PDF so the two match.
#
# Requires `typst` (>= 0.14) and `pdfinfo`/`pdftocairo` (poppler) on the PATH.

fit_journal_pdf <- function(pdf, max_w_mm = 180, max_h_mm = 225,
                            png = NULL, dpi = 300) {
  page_mm <- function(f) {
    info <- system2("pdfinfo", shQuote(f), stdout = TRUE)
    pts <- regmatches(info, regexpr("[0-9.]+ x [0-9.]+", info))[1]
    as.numeric(strsplit(pts, " x ")[[1]]) * 25.4 / 72
  }

  size <- page_mm(pdf)
  f <- min(1, max_w_mm / size[1], max_h_mm / size[2])

  if (f < 1) {
    tmp <- tempfile("fit_journal_")
    dir.create(tmp)
    file.copy(pdf, file.path(tmp, "in.pdf"))
    writeLines(c(
      sprintf("#set page(width: %.3fmm, height: %.3fmm, margin: 0pt)",
              size[1] * f, size[2] * f),
      '#image("in.pdf", width: 100%)'
    ), file.path(tmp, "fit.typ"))
    status <- system2("typst", c("compile", shQuote(file.path(tmp, "fit.typ")), shQuote(pdf)))
    unlink(tmp, recursive = TRUE)
    if (status != 0) stop("typst failed to fit ", pdf)
  }

  if (!is.null(png)) {
    system2("pdftocairo", c("-png", "-r", dpi, "-singlefile",
                            shQuote(pdf), shQuote(sub("\\.png$", "", png))))
  }

  new <- page_mm(pdf)
  message(sprintf("%s: %.0f x %.0f mm -> %.0f x %.0f mm",
                  basename(pdf), size[1], size[2], new[1], new[2]))
  invisible(new)
}
