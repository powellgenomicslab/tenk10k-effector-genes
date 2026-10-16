// Figure 6: Crohn's disease case study, formatted for Nature (180 mm wide, 8 pt bold
// panel labels), following scripts/1-overview/fig1.typ.
// Compile from the repo root (panels a-c are written by scripts/8-crohns/plot_figures/):
// typst compile --root . scripts/8-crohns/fig6.typ figures/fig-m6.pdf
//
// Layout: a, c and d stacked on the left; b, the tall heatmap, on the right.

#set page(width: 180mm, height: auto, margin: 2mm)
#set text(font: "Nimbus Sans", size: 8pt, weight: "bold")

#let dir = "/figures/panels/8-crohns/"

#let panel(tag, body) = block(width: 100%)[
  #body
  #place(top + left, dx: -1mm, tag)
]

// Panel b is drawn with fixed-aspect tiles on a 5 x 16 in page, so the heatmap
// fills only part of the page height. Show that part only: it starts 14.3% down
// the page and spans 71% of it (measured on the panel as written by
// 3-annotated_heatmap.R; re-measure if the number of genes or cell types changes).
#let b_top = 0.143
#let b_span = 0.71
#let b_ratio = 5 / 16

// Height-to-width ratios of the left-hand panels as saved (a: 5 x 2 in; c: 7 x 7 in;
// d: 5.5 x 5 in). The column widths are chosen so that the left column (a, c and d
// with 2 mm gaps) is exactly as tall as the visible part of b.
#let a_h = 2 / 5
#let c_h = 1
#let d_h = 5 / 5.5
#let gaps = 4mm
#let gutter = 3mm
#let usable = 180mm - 4mm
#let b_h = b_span / b_ratio
#let left_width = (b_h * (usable - gutter) - gaps) / (a_h + c_h + d_h + b_h)
#let right_width = usable - gutter - left_width
#let left_height = right_width * b_h
#let b_height = right_width / b_ratio

#grid(
  columns: (left_width, right_width),
  column-gutter: gutter,
  block(width: left_width, height: left_height, stack(
    spacing: 2mm,
    panel("a", image(dir + "fig6a_ext_source_annotation_barplot.pdf", width: 100%)),
    panel("c", image(dir + "fig6c_mr_deg_annot.pdf", width: 100%)),
    panel("d", image(dir + "fig6d_ibdverse_blood_enrichment.pdf", width: 100%)),
  )),
  panel("b", box(width: right_width, height: left_height, clip: true,
    place(top + left, dy: -b_height * b_top,
      image(dir + "fig6b_heatmap_canonical_deg_mrsenscoloc.pdf", width: right_width, height: b_height)))),
)
