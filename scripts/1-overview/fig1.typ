// Figure 1: study design, formatted for Nature (180 mm wide, 8 pt bold panel labels;
// Nimbus Sans, the metric equivalent of Helvetica installed in the figures environment).
// Compile from the repo root (panels b, c are written by scripts/1-overview/study_design.R):
// typst compile --root . scripts/1-overview/fig1.typ figures/fig-m1.pdf

#set page(width: 180mm, height: auto, margin: 2mm)
#set text(font: "Nimbus Sans", size: 8pt, weight: "bold")

#let panel(tag, path) = block(width: 100%)[
  #image(path, width: 100%)
  #place(top + left, dx: -1mm, tag)
]

#stack(spacing: 2mm,
  panel("a", "/figures/biorender/infographics_v5.pdf"),
  panel("b", "/figures/panels/1-overview/panel_b.svg"),
  panel("c", "/figures/panels/1-overview/panel_c.svg"),
)
