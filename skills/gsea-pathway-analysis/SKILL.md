---
name: gsea-pathway-analysis
description: Run pre-ranked GSEA on a differential-expression result table across one or more MSigDB collections (Hallmark, REACTOME, GO:BP, GO:MF, etc.) for a chosen species, and optionally drill into a specific gene-set signature with a running-enrichment plot and a per-sample ssGSEA activity heatmap. Use after a DESeq2 (or similar) DE analysis when the user wants pathway/gene-set enrichment, to test whether a known signature is enriched, or to see per-sample pathway activity. Step 2 of the DE → GSEA → report pipeline.
---

# GSEA pathway analysis

Pre-ranked GSEA with `clusterProfiler`, plus an optional signature deep-dive.

## When to use
- You have a DE results table and want **pathway enrichment** across MSigDB collections.
- You want to test whether a **specific signature** (a set of MSigDB gene-set IDs)
  is still enriched, and visualize its **per-sample activity** (ssGSEA).

## Inputs
- `<project>_DE_results.csv` from the DE skill (needs `GeneName` + a ranking column,
  default `stat`).
- For the signature ssGSEA heatmap: the original count matrix + metadata (same
  `config.R`), since per-sample scores need expression.

## How to run
1. Ensure `config.R` is present (shared with the DE skill) and `gene_set_collections`,
   `species`, `rank_metric` are set.
2. Copy `scripts/run_gsea.R` (and `scripts/explore_signature.R` if drilling into a
   signature) into the project.
3. `Rscript run_gsea.R`
4. Optional: set `signature_terms` (+ `signature_collection`/`subcollection`) in
   `config.R`, then `Rscript explore_signature.R`.

## Outputs (in `output_dir`)
Per collection `<NAME>`:
- `<project>_GSEA_<NAME>.csv`, `_dotplot.(pdf|png)`, `_barplot.(pdf|png)`.

Signature deep-dive:
- `<project>_<sig>_signature.csv` — NES/p/padj + ranks for the chosen gene sets.
- `<project>_<sig>_gseaplot.(pdf|png)` — running-enrichment plot with an NES/p table.
- `<project>_<sig>_ssgsea_scores.csv` and `_ssgsea_heatmap.(pdf|png)` — per-sample activity.

## Key conventions
- **Ranking metric** = DESeq2 `stat` by default; duplicate symbols kept by max |stat|.
- **Direction:** positive NES = enriched in the DE **test** group (matches the DE
  contrast direction).
- Compact, panel-ready figures: wrapped pathway labels, larger fonts (see the
  `figures-and-report` skill).
- Hallmark (`H`) has no subcollection — leave `subcollection = NA`.

## Requirements
R packages: `clusterProfiler`, `msigdbr`, `enrichplot`, `DOSE`, an org DB
(e.g. `org.Hs.eg.db`), `GSVA` (for ssGSEA), `dplyr`, `ggplot2`, plus `DESeq2`/`pheatmap`
for the signature heatmap.

Next step: **`figures-and-report`** to compile a PDF.
