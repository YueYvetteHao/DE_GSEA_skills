# =============================================================================
# run_gsea.R  —  Pre-ranked GSEA over multiple MSigDB collections
# Author: Yue Hao
# -----------------------------------------------------------------------------
# Ranks genes from a DE results table (config$rank_metric) and runs GSEA against
# each collection in config$gene_set_collections. Positive NES = enriched in the
# DE test group. Compact, panel-ready dotplot + bar plot per collection.
#
# Run (after the DE skill):  Rscript run_gsea.R
# =============================================================================

source("config.R")
cfg <- config

suppressPackageStartupMessages({
  library(clusterProfiler); library(msigdbr); library(enrichplot)
  library(DOSE); library(dplyr); library(ggplot2)
})
dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Ranked gene list from the DE table --------------------------------------
de_csv <- file.path(cfg$output_dir, paste0(cfg$project_name, "_DE_results.csv"))
if (!file.exists(de_csv)) stop("DE table not found: ", de_csv, "\nRun run_deseq2.R first.")
de <- read.csv(de_csv, check.names = FALSE)
m <- cfg$rank_metric
de <- de[!is.na(de[[m]]) & !is.na(de$GeneName) & nzchar(de$GeneName), ]
de <- de[order(abs(de[[m]]), decreasing = TRUE), ]
de <- de[!duplicated(de$GeneName), ]
gene_list <- sort(setNames(de[[m]], de$GeneName), decreasing = TRUE)

# --- Shared helper: wrap long pathway labels onto multiple lines --------------
wrap_labels <- function(x, width = 38)
  vapply(x, function(s) paste(strwrap(gsub("_", " ", s), width = width), collapse = "\n"), character(1))

plot_gsea_bar <- function(g, title, showCategory = 20, split = ".sign") {
  d <- as.data.frame(g)
  if (!".sign" %in% colnames(d)) d$.sign <- ifelse(d$NES >= 0, "positive", "negative")
  n <- floor(showCategory / 2)
  pos <- d[d$.sign == "positive", ]; neg <- d[d$.sign == "negative", ]
  pos <- pos[order(abs(pos$NES), decreasing = TRUE)[seq_len(min(n, nrow(pos)))], ]
  neg <- neg[order(abs(neg$NES), decreasing = TRUE)[seq_len(min(n, nrow(neg)))], ]
  d <- rbind(pos, neg)
  d$Description <- stats::reorder(d$Description, d$NES)
  ggplot(d, aes(x = Description, y = NES, fill = p.adjust)) +
    geom_col(width = 0.8) + coord_flip() +
    scale_x_discrete(labels = wrap_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_gradient(low = "red", high = "blue", name = "Adjusted p") +
    labs(x = NULL, y = "Normalized enrichment score (NES)", title = title) +
    theme_bw(base_size = 13) +
    theme(axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank())
}

run_one <- function(coll) {
  label <- paste0(cfg$project_name, ": ", coll$name)
  message("\n=== GSEA: ", coll$name, " ===")
  args <- list(species = cfg$species, collection = coll$collection)
  if (!is.na(coll$subcollection) && nzchar(coll$subcollection)) args$subcollection <- coll$subcollection
  t2g <- do.call(msigdbr, args) %>% dplyr::select(gs_name, gene_symbol)

  set.seed(42)
  g <- GSEA(geneList = gene_list, TERM2GENE = t2g, pvalueCutoff = 1,
            pAdjustMethod = "fdr", eps = 0, seed = TRUE)
  if (is.null(g) || nrow(as.data.frame(g)) == 0) { message("  no enriched sets."); return(invisible()) }

  prefix <- file.path(cfg$output_dir, paste0(cfg$project_name, "_GSEA_", coll$name))
  write.csv(as.data.frame(g), paste0(prefix, ".csv"), row.names = FALSE)

  dp <- dotplot(g, showCategory = 8, split = ".sign", font.size = 11, label_format = 35) +
    facet_grid(. ~ .sign) + ggtitle(label) + theme(plot.title = element_text(face = "bold"))
  ggsave(paste0(prefix, "_dotplot.pdf"), dp, width = 9, height = 6.5)
  ggsave(paste0(prefix, "_dotplot.png"), dp, width = 9, height = 6.5, dpi = 300)

  bp <- plot_gsea_bar(g, title = label, showCategory = 20)
  ggsave(paste0(prefix, "_barplot.pdf"), bp, width = 8.5, height = 8)
  ggsave(paste0(prefix, "_barplot.png"), bp, width = 8.5, height = 8, dpi = 300)
  message("  wrote table + dotplot + bar plot (", nrow(as.data.frame(g)), " sets).")
}

for (coll in cfg$gene_set_collections) run_one(coll)
message("\nrun_gsea.R complete.")
