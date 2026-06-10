# =============================================================================
# run_deseq2.R  —  Generalized two-group DESeq2 differential expression
# Author: Yue Hao
# -----------------------------------------------------------------------------
# Reads a count matrix + sample metadata (config.R), defines two groups by a
# threshold on a numeric column OR by an existing category column, and runs
# DESeq2 (test vs reference). Positive log2FC = higher in the test group.
#
# Run:  Rscript run_deseq2.R
# =============================================================================

source("config.R")
cfg <- config

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(readxl)
})
dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read any table type (.xlsx/.csv/tsv), stripping a BOM from headers -------
read_any <- function(f) {
  ext <- tolower(tools::file_ext(f))
  d <- if (ext %in% c("xlsx", "xls")) as.data.frame(readxl::read_excel(f), check.names = FALSE)
       else if (ext == "csv")         read.csv(f, check.names = FALSE)
       else                           read.delim(f, check.names = FALSE)
  colnames(d) <- sub("^﻿", "", colnames(d))
  d
}

# -----------------------------------------------------------------------------
# 1. Sample metadata -> two groups
# -----------------------------------------------------------------------------
meta <- read_any(cfg$metadata_file)
meta$.sample <- as.character(meta[[cfg$sample_id_col]])

if (cfg$group_mode == "threshold") {
  vcol <- grep(cfg$value_col_pattern, colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(vcol)) stop("No value column matching '", cfg$value_col_pattern, "'.")
  meta$.value <- suppressWarnings(as.numeric(meta[[vcol]]))
  meta <- meta[!is.na(meta$.value) & nzchar(meta$.sample), ]
  meta$.group <- ifelse(meta$.value <= cfg$threshold, cfg$low_label, cfg$high_label)
} else if (cfg$group_mode == "category") {
  meta$.group <- as.character(meta[[cfg$group_col]])
  meta <- meta[!is.na(meta$.group) & nzchar(meta$.sample), ]
  if (length(cfg$group_levels) == 2) meta <- meta[meta$.group %in% cfg$group_levels, ]
} else stop("group_mode must be 'threshold' or 'category'.")

groups <- unique(meta$.group)
if (length(groups) != 2) stop("Expected exactly 2 groups, got: ", paste(groups, collapse = ", "))
ref  <- cfg$reference_group
test <- setdiff(groups, ref)
if (length(test) != 1) stop("reference_group '", ref, "' is not one of the two groups.")

# -----------------------------------------------------------------------------
# 2. Count matrix
# -----------------------------------------------------------------------------
counts_raw  <- read.delim(cfg$count_file, check.names = FALSE)
gene_meta   <- counts_raw[, intersect(cfg$meta_cols, colnames(counts_raw)), drop = FALSE]
sample_cols <- setdiff(colnames(counts_raw), cfg$meta_cols)

# -----------------------------------------------------------------------------
# 3. Match samples to count columns
# -----------------------------------------------------------------------------
key <- function(x) if (cfg$match_by == "numeric_id") gsub("\\D", "", x) else x
expr_key <- key(sample_cols)
meta_key <- key(meta$.sample)
shared   <- intersect(expr_key, meta_key)
if (!length(shared)) stop("No samples shared between counts and metadata (check match_by).")

keep_cols <- sample_cols[expr_key %in% shared]
ord <- if (cfg$match_by == "numeric_id") order(as.numeric(gsub("\\D", "", keep_cols))) else order(keep_cols)
keep_cols <- keep_cols[ord]
meta_keep <- meta[match(key(keep_cols), meta_key), ]

# -----------------------------------------------------------------------------
# 4. Build matrix, pre-filter
# -----------------------------------------------------------------------------
cts <- as.matrix(counts_raw[, keep_cols, drop = FALSE])
storage.mode(cts) <- "integer"
rownames(cts) <- gene_meta[[cfg$gene_id_col]]
id2symbol <- setNames(gene_meta[[cfg$gene_symbol_col]], gene_meta[[cfg$gene_id_col]])
cts <- cts[rowSums(cts) >= cfg$min_count, , drop = FALSE]

# -----------------------------------------------------------------------------
# 5. DESeq2 (test vs reference)
# -----------------------------------------------------------------------------
coldata <- data.frame(Group = factor(meta_keep$.group, levels = c(ref, test)),
                      row.names = keep_cols)
message(sprintf("Samples matched: %d  (%s: %d, %s: %d)  | contrast: %s vs %s",
                nrow(coldata), test, sum(coldata$Group == test),
                ref, sum(coldata$Group == ref), test, ref))
if (any(table(coldata$Group) < 2)) stop("Need >= 2 samples per group.")

dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ Group)
dds$Group <- relevel(dds$Group, ref = ref)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", test, ref), alpha = 0.05)

res_df <- as.data.frame(res)
res_df$GeneId   <- rownames(res_df)
res_df$GeneName <- id2symbol[res_df$GeneId]
res_df <- res_df[, c("GeneName", "GeneId", "baseMean", "log2FoldChange",
                     "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]

out_csv <- file.path(cfg$output_dir, paste0(cfg$project_name, "_DE_results.csv"))
write.csv(res_df, out_csv, row.names = FALSE)
message("Wrote DE table: ", out_csv, "  (positive log2FC = higher in ", test, ")")
print(summary(res))

# -----------------------------------------------------------------------------
# 6. Volcano (gene symbols)
# -----------------------------------------------------------------------------
v <- EnhancedVolcano(res_df, lab = res_df$GeneName, x = "log2FoldChange", y = "pvalue",
                     title = paste0(cfg$project_name, ": ", test, " vs ", ref),
                     subtitle = paste0("Positive log2FC = higher in ", test), labSize = 3.5)
ggplot2::ggsave(file.path(cfg$output_dir, paste0(cfg$project_name, "_volcano.pdf")), v, width = 9, height = 8)
ggplot2::ggsave(file.path(cfg$output_dir, paste0(cfg$project_name, "_volcano.png")), v, width = 9, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 7. Top-40 heatmap (20 up + 20 down), VST, z-scored, symbol rows
# -----------------------------------------------------------------------------
ranked <- res_df[!is.na(res_df$padj), ]
ranked <- ranked[order(ranked$log2FoldChange, decreasing = TRUE), ]
top_genes <- rbind(head(ranked, 20), tail(ranked, 20))

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes$GeneId, , drop = FALSE]
rownames(mat) <- top_genes$GeneName
if (cfg$match_by == "numeric_id") colnames(mat) <- meta_keep$.sample  # show sample names

# Column annotation: Group (+ z-scored value in threshold mode).
annotation_col <- data.frame(Group = coldata$Group, row.names = colnames(mat))
ann_colors <- list(Group = setNames(c("#F8766D", "#619CFF"), c(test, ref)))
if (cfg$group_mode == "threshold") {
  v_raw <- meta_keep$.value
  z <- as.numeric(scale(if (all(v_raw > 0)) log10(v_raw) else v_raw))
  annotation_col$value_z <- z
  ann_colors$value_z <- c("#FFFFD9", "#225EA8")
  cap <- paste0(cfg$low_label, ": ", cfg$value_col_pattern, " <= ", cfg$threshold,
                " ", cfg$value_units, "   |   ", cfg$high_label, ": > ",
                cfg$threshold, " ", cfg$value_units)
} else {
  cap <- paste0("Groups: ", test, " vs ", ref, " (reference)")
}

draw <- function(file) {
  pheatmap(mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_col = annotation_col, annotation_colors = ann_colors,
           border_color = NA, color = colorRampPalette(c("navy", "white", "firebrick"))(50),
           cellwidth = 16, cellheight = 11, fontsize_row = 9, fontsize_col = 9, angle_col = 45,
           main = paste0(cfg$project_name, ": top 40 DE genes (20 up / 20 down)\n", cap),
           filename = file)
}
draw(file.path(cfg$output_dir, paste0(cfg$project_name, "_top40_heatmap.pdf")))
draw(file.path(cfg$output_dir, paste0(cfg$project_name, "_top40_heatmap.png")))
message("run_deseq2.R complete.")
