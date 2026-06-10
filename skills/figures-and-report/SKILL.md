---
name: figures-and-report
description: Apply consistent, publication/panel-ready figure styling to RNA-seq analysis plots and compile a PDF summary report (R Markdown) covering the data, sample/group table, DE and GSEA results, any signature analysis, and all key figures. Use when the user wants to optimize/standardize figures (compact, readable, larger fonts, wrapped labels, gene-symbol labels, group + value annotations on heatmaps) or to generate a shareable PDF/HTML report from existing analysis outputs. Step 3 (final) of the DE → GSEA → report pipeline.
---

# Figures & report

Centralizes figure conventions and builds a PDF report from the analysis outputs.

## When to use
- The user asks to **make figures more compact / readable / panel-ready**, fix
  small text, wrap long pathway names, or standardize plot style.
- The user wants a **PDF (or HTML) report** summarizing the analysis.

## Figure conventions (also embedded in the DE/GSEA scripts)
Provided in `scripts/figure_style.R` (source it for ad-hoc plots):
- `wrap_labels()` — wrap long `UNDERSCORE_NAMES` onto multiple lines.
- `theme_report()` — `theme_bw` with larger base font, bold title.
- `gsea_bar()` / heatmap defaults — compact, larger fonts, auto-sized.
- Rules: **gene symbols** as labels (never raw IDs); state the **contrast
  direction** in titles/captions; annotate heatmaps with the group (and value).
- Sizing: avoid over-wide canvases (they shrink text in panels); let `pheatmap`
  auto-size when using fixed `cellwidth`/`cellheight`; widen only enough to fit
  wrapped labels.

## Report
`scripts/report_template.Rmd` + `scripts/render_report.R` produce
`outputs/<project>_report.pdf`. The report reads `config.R` and the files written
by the DE/GSEA skills, so it adapts to the project automatically (threshold vs
category mode, whichever collections were run, optional signature section).

## How to run
1. Copy `scripts/report_template.Rmd`, `scripts/render_report.R`, and
   `scripts/figure_style.R` into the project.
2. `Rscript render_report.R`

## Requirements
- R: `rmarkdown`, `knitr`, `readxl`.
- **pandoc** (bundled with RStudio; `render_report.R` auto-detects it on macOS) and a
  **LaTeX** engine for PDF: `Rscript -e 'tinytex::install_tinytex()'` (one-time).
- In RStudio you can also open `report_template.Rmd` and click **Knit**.
