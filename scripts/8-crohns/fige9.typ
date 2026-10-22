// Extended Data Figure 9: GPX1 at the Crohn's disease locus, 180 x 195 mm as published,
// with 8 pt bold panel labels as in fig6.typ.
// Compile from the repo root (panels are written by plot_figures/2-plot_locus_zoom_fig.R):
// typst compile --root . scripts/8-crohns/fige9.typ figures/fig-e9.pdf
//
// Layout: a, the three locus zoom plots, on the left; b and c, the cross plots, on the right.

#set page(width: 180mm, height: 195mm, margin: 2mm)
#set text(font: "Nimbus Sans", size: 8pt, weight: "bold")

#let dir = "/figures/panels/8-crohns/"

#let panel(tag, body) = block[
  #body
  #place(top + left, tag)
]

#grid(
  columns: (88mm, 88mm),
  panel("a", image(dir + "fige9a_locus_zoom.pdf", width: 88mm, height: 191mm)),
  stack(
    spacing: 3mm,
    panel("b", image(dir + "fige9b_cross_plot.pdf", width: 88mm, height: 94mm)),
    panel("c", image(dir + "fige9c_cross_plot.pdf", width: 88mm, height: 94mm)),
  ),
)
