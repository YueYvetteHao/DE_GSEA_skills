---
name: deseq2-differential-expression
description: Run a two-group RNA-seq differential expression analysis with DESeq2 from a raw count matrix and a sample metadata table. Use when the user wants to compare gene expression between two groups of samples/cell lines (e.g. sensitive vs resistant, treated vs control, mutant vs WT) — including when groups are defined by thresholding a numeric column (IC50, dose) or by an existing category column. Produces a DESeq2 results table (gene symbols), a volcano plot, and a top-N differentially expressed gene heatmap. Step 1 of the DE → GSEA → report pipeline.
---

# DESeq2 two-group differential expression

Generalized DESeq2 workflow. Compares a **test** group vs a **reference** group and
writes a symbol-labeled results table plus a volcano and a top-N heatmap.

## When to use
- Any two-group bulk RNA-seq comparison from a raw **count matrix**.
- Groups defined either by a **threshold** on a numeric column (e.g. IC50 ≤ 10 →
  Sensitive) or by an **existing category** column (e.g. Treatment = drug/vehicle).

## Inputs
1. **Count matrix** (tab-delimited): some gene-metadata columns (feature id, symbol,
   length, …) followed by one column per sample.
2. **Sample metadata** (`.xlsx`/`.csv`): one row per sample, with a sample-name
   column and either a numeric column (threshold mode) or a group column (category mode).

Samples are matched to count columns by **numeric ID** (e.g. `X14` ≡ `Sample14`) or by
**exact name** — set `match_by` in the config.

## How to run
1. Copy `scripts/config.R` and `scripts/run_deseq2.R` into the project directory.
2. Edit `config.R` — at minimum `count_file`, `metadata_file`, the group definition
   (`group_mode` + its fields), `reference_group`, and `project_name`.
3. `Rscript run_deseq2.R`

## Outputs (in `output_dir`)
- `<project>_DE_results.csv` — full DESeq2 table, **GeneName first**, then feature id and stats.
- `<project>_volcano.(pdf|png)` — labeled by gene symbol.
- `<project>_top40_heatmap.(pdf|png)` — top up/down genes, columns annotated by group
  (and z-scored value in threshold mode), with the group-definition caption.

## Key conventions
- **Direction:** contrast = `test vs reference` (reference = baseline). Positive
  log2FC ⇒ higher in the test group. In threshold mode the low-value group is the
  test group by default (set `reference_group` to the high-value group).
- **Symbols everywhere;** feature IDs are internal/secondary.
- Genes with total count `< min_count` are pre-filtered.

## Requirements
R packages: `DESeq2`, `EnhancedVolcano`, `pheatmap`, `readxl` (Bioconductor for
DESeq2). See the `figures-and-report` skill for figure conventions reused here.

Next step: **`gsea-pathway-analysis`** (reads `<project>_DE_results.csv`).
