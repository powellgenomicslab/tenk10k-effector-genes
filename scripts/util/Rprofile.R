# Read by every R step of scripts/Snakefile (via R_PROFILE_USER).

# A script that prints a plot at the top level makes R open its default device,
# which writes Rplots.pdf into the working directory; draw such plots on a null
# PDF device instead (same font metrics, no file). Figures are saved with ggsave
# as before.
options(device = function(...) grDevices::pdf(NULL))

# Text that names no font uses the generic family "sans", which fontconfig maps to
# DejaVu Sans in the figures environment. Prefer Helvetica, as in the manuscript,
# and Nimbus Sans (a metric-compatible Helvetica clone installed in the environment)
# where Helvetica is not installed. The rule goes in a per-user fontconfig file that
# the environment's configuration reads (conf.d/50-user.conf), so it applies to
# ragg (PNG), cairo_pdf and pdftoppm alike; it must be in place before any of them
# first loads fonts. On macOS, ragg finds Helvetica without fontconfig.
local({
  home <- file.path(tempdir(), "xdg")
  dir.create(file.path(home, "fontconfig"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    '<?xml version="1.0"?>',
    '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">',
    "<fontconfig>",
    sprintf(paste0('  <alias binding="strong"><family>%s</family>',
                   "<prefer><family>Helvetica</family><family>Nimbus Sans</family></prefer></alias>"),
            c("sans-serif", "sans")),
    "</fontconfig>"
  ), file.path(home, "fontconfig", "fonts.conf"))
  Sys.setenv(XDG_CONFIG_HOME = home)
})
