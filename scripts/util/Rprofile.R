# Read by every R step of scripts/Snakefile (via R_PROFILE_USER). A script that
# prints a plot at the top level makes R open its default device, which writes
# Rplots.pdf into the working directory; draw such plots on a null PDF device
# instead (same font metrics, no file). Figures are saved with ggsave as before.
options(device = function(...) grDevices::pdf(NULL))
