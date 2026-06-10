# RNA-seq DE + GSEA Skills

A set of three composable, reusable **skills** for bulk RNA-seq analysis:
differential expression → gene-set enrichment → a publication-ready PDF report.
Each skill bundles generalized R scripts plus a `SKILL.md` that lets an AI coding
assistant (e.g. Claude Code) discover and run it — but the scripts are equally
usable on their own from the command line.

> Built for two-group comparisons (e.g. treated vs control, mutant vs wild-type,
> sensitive vs resistant). Groups can be defined by an existing category column or
> by thresholding a numeric column (dose, IC50, score).

## Skills

| Skill | Purpose |
|-------|---------|
| [`deseq2-differential-expression`](deseq2-differential-expression/) | Two-group differential expression with **DESeq2** → results table, volcano, top-N heatmap. |
| [`gsea-pathway-analysis`](gsea-pathway-analysis/) | Pre-ranked **GSEA** over multiple MSigDB collections + optional single-signature deep-dive (running-enrichment plot + per-sample ssGSEA heatmap). |
| [`figures-and-report`](figures-and-report/) | Shared figure-styling conventions + an adaptive **R Markdown PDF report**. |

## Features

- **Generalized** — any two-group comparison; groups by category *or* numeric threshold.
- **Flexible sample matching** — by numeric ID (e.g. `X14` ≡ `Sample14`) or exact name.
- **Configurable** — species, GSEA ranking metric, and MSigDB collections via one `config.R`.
- **Consistent conventions** — gene symbols everywhere, explicit contrast direction,
  compact panel-ready figures.
- **Reproducible** — one config file drives the whole pipeline; re-run to regenerate everything.

## Requirements

- **R ≥ 4.2**
- R packages:
  - Bioconductor: `DESeq2`, `clusterProfiler`, `enrichplot`, `DOSE`, an organism
    annotation DB (e.g. `org.Hs.eg.db`), `GSVA`
  - CRAN: `msigdbr`, `dplyr`, `ggplot2`, `EnhancedVolcano`, `pheatmap`, `readxl`,
    `rmarkdown`, `knitr`
- For the PDF report: **pandoc** (bundled with RStudio; auto-detected) and a
  **LaTeX** engine — install once with `Rscript -e 'tinytex::install_tinytex()'`.

Install the R packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","clusterProfiler","enrichplot","DOSE","org.Hs.eg.db","GSVA"))
install.packages(c("msigdbr","dplyr","ggplot2","EnhancedVolcano","pheatmap",
                   "readxl","rmarkdown","knitr"))
```

## How to use

### 1. Set up a project directory

Copy the scripts you need into a working directory alongside your data, and copy
the master config:

```sh
mkdir my_project && cd my_project
cp /path/to/skills/deseq2-differential-expression/scripts/config.R .
cp /path/to/skills/deseq2-differential-expression/scripts/run_deseq2.R .
cp /path/to/skills/gsea-pathway-analysis/scripts/run_gsea.R .
cp /path/to/skills/gsea-pathway-analysis/scripts/explore_signature.R .   # optional
cp /path/to/skills/figures-and-report/scripts/{report_template.Rmd,render_report.R,figure_style.R} .
```

### 2. Prepare two inputs

- **Count matrix** (tab-delimited): gene-metadata columns (feature id, symbol,
  length …) followed by one column per sample.

  ```
  GeneId           GeneName  Length  X1    X2    X3   ...
  ENSG00000123456  TP53      2500    1043  876   1290 ...
  ```

- **Sample metadata** (`.xlsx` or `.csv`): one row per sample.

  ```
  Sample,Treatment,Dose
  Sample1,Drug,12.5
  Sample2,Vehicle,0.4
  ```

### 3. Edit `config.R`

Set the inputs, the group definition, and a project name. Two ways to define groups:

```r
# (a) THRESHOLD mode — derive groups from a numeric column
group_mode        = "threshold",
value_col_pattern = "Dose",      # regex to find the numeric column
threshold         = 10,          # value <= threshold -> low_label
low_label         = "Low",
high_label        = "High",
reference_group   = "High",      # baseline; positive log2FC = higher in "Low"

# (b) CATEGORY mode — use an existing two-level column
group_mode      = "category",
group_col       = "Treatment",
group_levels    = c("Drug","Vehicle"),
reference_group = "Vehicle",     # positive log2FC = higher in "Drug"
```

### 4. Run the pipeline

```sh
Rscript run_deseq2.R          # -> <project>_DE_results.csv, volcano, top-40 heatmap
Rscript run_gsea.R            # -> per-collection GSEA tables + dotplots/barplots
Rscript explore_signature.R   # optional: set config$signature_terms first
Rscript render_report.R       # -> <project>_report.pdf
```

All outputs land in `config$output_dir` (default `outputs/`). In RStudio you can
also open `report_template.Rmd` and click **Knit**.

### Using with an AI assistant (Claude Code)

Copy a skill folder into `~/.claude/skills/` (user) or `<project>/.claude/skills/`
(project). The assistant reads each `SKILL.md` and invokes the skill when your
request matches its description (e.g. “run DESeq2 on these counts”, “run GSEA”,
“make a PDF report”).

## Outputs

| File | From | Description |
|------|------|-------------|
| `<project>_DE_results.csv` | DE | DESeq2 table, gene symbol first |
| `<project>_volcano.(pdf\|png)` | DE | Volcano plot (symbol labels) |
| `<project>_top40_heatmap.(pdf\|png)` | DE | Top 20 up / 20 down genes, annotated |
| `<project>_GSEA_<NAME>.csv` + `_dotplot`/`_barplot` | GSEA | Per-collection results & figures |
| `<project>_<sig>_signature.csv`, `_gseaplot`, `_ssgsea_*` | GSEA | Optional signature deep-dive |
| `<project>_report.pdf` | report | Full summary report |

## Conventions (enforced across all skills)

- **Direction is fixed and explicit:** contrast = `test vs reference`; positive
  log2FC / NES = higher in the test group. Stated in tables, captions, and the report.
- **Gene symbols** label every table and figure; feature IDs stay as a secondary column.
- **GSEA ranking:** DESeq2 `stat` by default; duplicate symbols resolved by max |stat|.
- **Compact, panel-ready figures:** wrapped pathway labels, larger fonts, auto-sized heatmaps.

## Repository layout

```
skills/
├── README.md
├── deseq2-differential-expression/
│   ├── SKILL.md
│   └── scripts/{config.R, run_deseq2.R}
├── gsea-pathway-analysis/
│   ├── SKILL.md
│   └── scripts/{run_gsea.R, explore_signature.R}
└── figures-and-report/
    ├── SKILL.md
    └── scripts/{figure_style.R, report_template.Rmd, render_report.R}
```

## Author & license

Created by **Yue Hao**. Released for public use — add a `LICENSE` file (e.g. MIT)
before publishing if you want to set explicit terms.

Contributions and issues welcome.
