# =============================================================================
# figure_style.R  —  Shared figure conventions for the DE/GSEA/report skills
# Author: Yue Hao
# -----------------------------------------------------------------------------
# Source this for consistent, panel-ready plots:  source("figure_style.R")
# Conventions: gene SYMBOLS as labels (never raw IDs); state contrast direction
# in titles; wrap long pathway names; larger fonts; auto-sized heatmaps.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# Wrap long UNDERSCORE_NAMES (REACTOME_/GOBP_/...) onto multiple lines so labels
# stay legible without a very wide canvas.
wrap_labels <- function(x, width = 38) {
  vapply(x, function(s) paste(strwrap(gsub("_", " ", s), width = width),
                              collapse = "\n"), character(1))
}

# Compact base theme: larger fonts, bold title, light grid.
theme_report <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = max(9, base_size - 3)),
          panel.grid.minor = element_blank())
}

# NES bar plot from a GSEA result (top |NES| per direction). Compact + wrapped.
gsea_bar <- function(gsea_result, title, showCategory = 20) {
  d <- as.data.frame(gsea_result)
  if (!".sign" %in% colnames(d)) d$.sign <- ifelse(d$NES >= 0, "positive", "negative")
  n <- floor(showCategory / 2)
  pos <- d[d$.sign == "positive", ]; neg <- d[d$.sign == "negative", ]
  pos <- pos[order(abs(pos$NES), decreasing = TRUE)[seq_len(min(n, nrow(pos)))], ]
  neg <- neg[order(abs(neg$NES), decreasing = TRUE)[seq_len(min(n, nrow(neg)))], ]
  d <- rbind(pos, neg); d$Description <- stats::reorder(d$Description, d$NES)
  ggplot(d, aes(x = Description, y = NES, fill = p.adjust)) +
    geom_col(width = 0.8) + coord_flip() +
    scale_x_discrete(labels = wrap_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_gradient(low = "red", high = "blue", name = "Adjusted p") +
    labs(x = NULL, y = "Normalized enrichment score (NES)", title = title) +
    theme_report() + theme(panel.grid.major.y = element_blank())
}

# Recommended pheatmap arguments (pass via do.call or copy inline). Fixed cell
# size + no width/height lets pheatmap auto-size so labels are never clipped.
heatmap_defaults <- list(
  scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

# Save a ggplot to both PDF and PNG with one call.
save_fig <- function(plot, path_noext, width = 9, height = 7) {
  ggsave(paste0(path_noext, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(path_noext, ".png"), plot, width = width, height = height, dpi = 300)
}
