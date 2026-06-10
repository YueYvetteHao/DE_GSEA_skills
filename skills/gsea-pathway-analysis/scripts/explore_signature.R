# =============================================================================
# explore_signature.R  —  Drill into a specific gene-set signature
# Author: Yue Hao
# -----------------------------------------------------------------------------
# For a user-specified set of MSigDB gene-set IDs (config$signature_terms):
#   A. group-level GSEA significance + a running-enrichment (gseaplot) figure
#   B. per-sample ssGSEA activity heatmap (needs the count matrix + metadata)
#
# Run (after run_deseq2.R):  Rscript explore_signature.R
# =============================================================================

source("config.R")
cfg <- config
if (length(cfg$signature_terms) == 0) stop("Set config$signature_terms to one or more MSigDB gene-set IDs.")

suppressPackageStartupMessages({
  library(DESeq2); library(clusterProfiler); library(msigdbr)
  library(enrichplot); library(GSVA); library(pheatmap); library(dplyr); library(ggplot2); library(readxl)
})
dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)
sig <- gsub("[^A-Za-z0-9]+", "_", cfg$signature_label)
terms <- cfg$signature_terms

read_any <- function(f) {
  ext <- tolower(tools::file_ext(f))
  d <- if (ext %in% c("xlsx","xls")) as.data.frame(readxl::read_excel(f), check.names = FALSE)
       else if (ext == "csv") read.csv(f, check.names = FALSE) else read.delim(f, check.names = FALSE)
  colnames(d) <- sub("^﻿", "", colnames(d)); d
}

# =============================================================================
# PART A — group-level GSEA on the signature's collection
# =============================================================================
de_csv <- file.path(cfg$output_dir, paste0(cfg$project_name, "_DE_results.csv"))
if (!file.exists(de_csv)) stop("DE table not found: ", de_csv)
de <- read.csv(de_csv, check.names = FALSE)
m <- cfg$rank_metric
de <- de[!is.na(de[[m]]) & !is.na(de$GeneName) & nzchar(de$GeneName), ]
de <- de[order(abs(de[[m]]), decreasing = TRUE), ]; de <- de[!duplicated(de$GeneName), ]
gene_list <- sort(setNames(de[[m]], de$GeneName), decreasing = TRUE)

args <- list(species = cfg$species, collection = cfg$signature_collection)
if (!is.na(cfg$signature_subcollection) && nzchar(cfg$signature_subcollection))
  args$subcollection <- cfg$signature_subcollection
t2g <- do.call(msigdbr, args) %>% dplyr::select(gs_name, gene_symbol)

set.seed(42)
gsea <- GSEA(gene_list, TERM2GENE = t2g, pvalueCutoff = 1, pAdjustMethod = "fdr", eps = 0, seed = TRUE)
res <- as.data.frame(gsea)
res$rank_by_padj <- rank(res$p.adjust, ties.method = "min")
present <- terms[terms %in% res$ID]
if (!length(present)) stop("None of signature_terms found in the GSEA output.")

retro <- res[res$ID %in% present, c("ID", "setSize", "NES", "pvalue", "p.adjust", "rank_by_padj")]
retro <- retro[order(retro$p.adjust), ]
cat("\n==== ", cfg$signature_label, " ====\n"); print(retro, row.names = FALSE)
cat("All significant (padj<0.05): ", all(retro$p.adjust < 0.05), "\n")
write.csv(retro, file.path(cfg$output_dir, paste0(cfg$project_name, "_", sig, "_signature.csv")), row.names = FALSE)

# Running-enrichment plot with NES + p table (NES added via enrichplot helper).
ep <- gseaplot2(gsea, geneSetID = present,
                title = paste0(cfg$project_name, ": ", cfg$signature_label), pvalue_table = FALSE)
pd <- res[match(present, res$ID), c("Description", "NES", "pvalue", "p.adjust")]
rownames(pd) <- pd$Description; pd <- pd[, -1, drop = FALSE]
pd$NES <- format(round(pd$NES, 2), nsmall = 2); pd$pvalue <- format(pd$pvalue, digits = 3)
pd$p.adjust <- format(pd$p.adjust, digits = 3)
tp <- enrichplot:::tableGrob2(pd, ep[[1]])
ep[[1]] <- ep[[1]] + theme(legend.position = "none") + coord_cartesian(clip = "off") +
  annotation_custom(tp, xmin = quantile(ep[[1]]$data$x, 0.5), xmax = quantile(ep[[1]]$data$x, 0.95),
                    ymin = quantile(ep[[1]]$data$runningScore, 0.7), ymax = quantile(ep[[1]]$data$runningScore, 0.9))
for (i in seq_along(ep)) ep[[i]] <- ep[[i]] + theme(plot.margin = margin(5.5, 150, 5.5, 5.5))
ggsave(file.path(cfg$output_dir, paste0(cfg$project_name, "_", sig, "_gseaplot.pdf")), ep, width = 14, height = 8)
ggsave(file.path(cfg$output_dir, paste0(cfg$project_name, "_", sig, "_gseaplot.png")), ep, width = 14, height = 8, dpi = 300)

# =============================================================================
# PART B — per-sample ssGSEA activity heatmap
# =============================================================================
meta <- read_any(cfg$metadata_file); meta$.sample <- as.character(meta[[cfg$sample_id_col]])
if (cfg$group_mode == "threshold") {
  vcol <- grep(cfg$value_col_pattern, colnames(meta), ignore.case = TRUE, value = TRUE)[1]
  meta$.value <- suppressWarnings(as.numeric(meta[[vcol]]))
  meta <- meta[!is.na(meta$.value), ]
  meta$.group <- ifelse(meta$.value <= cfg$threshold, cfg$low_label, cfg$high_label)
} else { meta$.group <- as.character(meta[[cfg$group_col]]); meta <- meta[!is.na(meta$.group), ] }

counts_raw <- read.delim(cfg$count_file, check.names = FALSE)
gene_meta  <- counts_raw[, intersect(cfg$meta_cols, colnames(counts_raw)), drop = FALSE]
sample_cols <- setdiff(colnames(counts_raw), cfg$meta_cols)
key <- function(x) if (cfg$match_by == "numeric_id") gsub("\\D", "", x) else x
shared <- intersect(key(sample_cols), key(meta$.sample))
keep_cols <- sample_cols[key(sample_cols) %in% shared]
keep_cols <- keep_cols[if (cfg$match_by == "numeric_id") order(as.numeric(gsub("\\D","",keep_cols))) else order(keep_cols)]
meta_keep <- meta[match(key(keep_cols), key(meta$.sample)), ]

cts <- as.matrix(counts_raw[, keep_cols, drop = FALSE]); storage.mode(cts) <- "integer"
rownames(cts) <- gene_meta[[cfg$gene_id_col]]
cts <- cts[rowSums(cts) >= cfg$min_count, , drop = FALSE]
id2symbol <- setNames(gene_meta[[cfg$gene_symbol_col]], gene_meta[[cfg$gene_id_col]])

coldata <- data.frame(Group = factor(meta_keep$.group), row.names = keep_cols)
vsd <- vst(DESeqDataSetFromMatrix(cts, coldata, design = ~ 1), blind = TRUE)
expr <- assay(vsd); sym <- id2symbol[rownames(expr)]
ok <- !is.na(sym) & nzchar(sym); expr <- expr[ok, ]; sym <- sym[ok]
o <- order(rowMeans(expr), decreasing = TRUE); expr <- expr[o, ]; sym <- sym[o]
keep <- !duplicated(sym); expr <- expr[keep, ]; rownames(expr) <- sym[keep]
colnames(expr) <- meta_keep$.sample

gs <- t2g %>% dplyr::filter(gs_name %in% present)
gene_sets <- split(gs$gene_symbol, gs$gs_name)
scores <- gsva(ssgseaParam(exprData = expr, geneSets = gene_sets))
rownames(scores) <- sub(paste0("^", cfg$signature_collection, "_|^REACTOME_|^GOBP_|^GOMF_|^HALLMARK_"), "", rownames(scores))
write.csv(as.data.frame(scores), file.path(cfg$output_dir, paste0(cfg$project_name, "_", sig, "_ssgsea_scores.csv")))

annotation_col <- data.frame(Group = coldata$Group, row.names = colnames(expr))
draw <- function(file) pheatmap(scores, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = annotation_col, border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  cellwidth = 20, cellheight = 28, fontsize_row = 9, fontsize_col = 9, angle_col = 45,
  main = paste0(cfg$project_name, ": ", cfg$signature_label, " ssGSEA (per sample)"), filename = file)
draw(file.path(cfg$output_dir, paste0(cfg$project_name, "_", sig, "_ssgsea_heatmap.pdf")))
draw(file.path(cfg$output_dir, paste0(cfg$project_name, "_", sig, "_ssgsea_heatmap.png")))
message("explore_signature.R complete.")
