# =============================================================================
# config.R  —  Master settings for the DE + GSEA + report skills
# Author: Yue Hao
# -----------------------------------------------------------------------------
# This single file is read by all three skills. Edit it per project, then run
# run_deseq2.R -> run_gsea.R -> (explore_signature.R) -> render_report.R.
# =============================================================================

config <- list(

  # ---- Inputs -------------------------------------------------------------
  count_file    = "counts.txt",        # tab-delimited: gene-meta cols + sample cols
  metadata_file = "samples.xlsx",      # per-sample table (.xlsx/.xls/.csv)

  # ---- Sample <-> count-column matching -----------------------------------
  sample_id_col = "Models",            # metadata column holding sample names
  match_by      = "numeric_id",        # "numeric_id" (X14 == Sample14) or "exact"

  # ---- Group definition ---------------------------------------------------
  # Two ways to define the two groups; pick one with group_mode.
  group_mode    = "threshold",         # "threshold" or "category"

  # (a) threshold mode: derive groups from a numeric column
  #     (e.g. a dose, IC50, or any per-sample score)
  value_col_pattern = "Dose",          # regex to find the numeric column (case-insensitive)
  threshold     = 10,                  # value <= threshold -> low_label
  low_label     = "Low",               # group for value <= threshold
  high_label    = "High",              # group for value >  threshold
  value_units   = "units",             # for figure/report captions only

  # (b) category mode: use an existing two-level column
  group_col     = NA,                  # e.g. "Treatment"
  group_levels  = NA,                  # e.g. c("Drug","Vehicle")

  # The reference (baseline) group. Positive log2FC/NES = higher in the OTHER
  # (test) group. In threshold mode this is usually the high-value group.
  reference_group = "High",

  # ---- Count matrix layout ------------------------------------------------
  meta_cols       = c("GeneId", "GeneName", "Length"),  # non-sample columns
  gene_id_col     = "GeneId",          # unique feature id (used for rownames)
  gene_symbol_col = "GeneName",        # display label everywhere
  min_count       = 10,                # drop genes with total count < this

  # ---- Output -------------------------------------------------------------
  project_name = "MyProject",
  output_dir   = "outputs",

  # ---- GSEA / MSigDB ------------------------------------------------------
  species     = "Homo sapiens",
  rank_metric = "stat",                # DESeq2 column to rank by (stat | log2FoldChange)
  gene_set_collections = list(
    list(name = "HALLMARK", collection = "H",  subcollection = NA),
    list(name = "REACTOME", collection = "C2", subcollection = "CP:REACTOME"),
    list(name = "GO_BP",    collection = "C5", subcollection = "GO:BP"),
    list(name = "GO_MF",    collection = "C5", subcollection = "GO:MF")
  ),

  # ---- Optional signature deep-dive (explore_signature.R) ------------------
  signature_label         = "Signature of interest",
  signature_collection    = "C2",
  signature_subcollection = "CP:REACTOME",
  signature_terms = c()                # MSigDB gene-set IDs to drill into, e.g.
                                       # c("REACTOME_CHOLESTEROL_BIOSYNTHESIS")
)
