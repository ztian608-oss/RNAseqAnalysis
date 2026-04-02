scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- 1
  (x - m) / s
}

plot_deseq2_pca <- function(vsd, output_pdf, output_csv, intgroup) {
  pca_data <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"), 2)
  pca_data$PC1_var_percent <- percent_var[1]
  pca_data$PC2_var_percent <- percent_var[2]
  pca_data$PC1_PC2_total_percent <- sum(percent_var[1:2])

  p <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, color = .data[[intgroup]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "PCA (DESeq2 VST)",
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    )

  ggplot2::ggsave(output_pdf, plot = p, width = 5, height = 4)
  readr::write_csv(pca_data, output_csv)
  pca_data
}

run_volcano <- function(tbl, nm, padj_cutoff, log2_cutoff, outdir) {
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    stop("EnhancedVolcano package is required for volcano plot generation.")
  }
  tbl <- tbl[!is.na(tbl$pvalue) & !is.na(tbl$log2FoldChange), , drop = FALSE]
  if (nrow(tbl) == 0) return(NULL)
  out_file <- file.path(outdir, paste0("volcano_", nm, ".pdf"))
  labels <- ifelse(is.na(tbl$gene_name) | tbl$gene_name == "", tbl$gene_id, tbl$gene_name)
  tbl$volcano_group <- dplyr::case_when(
    !is.na(tbl$padj) & tbl$padj <= padj_cutoff & tbl$log2FoldChange >= log2_cutoff ~ "DEG Up",
    !is.na(tbl$padj) & tbl$padj <= padj_cutoff & tbl$log2FoldChange <= -log2_cutoff ~ "DEG Down",
    TRUE ~ "Not Significant"
  )
  keyvals <- c("DEG Up" = "#b2182b", "DEG Down" = "#2166ac", "Not Significant" = "grey70")[tbl$volcano_group]
  names(keyvals) <- tbl$volcano_group
  p <- EnhancedVolcano::EnhancedVolcano(
    tbl,
    lab = labels,
    x = "log2FoldChange",
    y = "pvalue",
    title = paste("Volcano plot:", nm),
    subtitle = paste0("x = log2FoldChange; y = -log10(p value); point classes use padj <= ", padj_cutoff,
      " and |log2FC| >= ", log2_cutoff),
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~ italic(P)),
    pCutoff = 1,
    FCcutoff = 0,
    labSize = 3,
    pointSize = 1.8,
    colAlpha = 0.9,
    colCustom = keyvals
  )
  ggplot2::ggsave(out_file, plot = p, width = 5, height = 5)
  out_file
}

run_ma_plot <- function(tbl, nm, padj_cutoff, log2_cutoff, outdir) {
  if (!all(c("baseMean", "log2FoldChange", "padj") %in% colnames(tbl))) {
    return(NULL)
  }
  x <- tbl[!is.na(tbl$baseMean) & !is.na(tbl$log2FoldChange), , drop = FALSE]
  if (nrow(x) == 0) return(NULL)

  x$significance <- dplyr::case_when(
    !is.na(x$padj) & x$padj <= padj_cutoff & x$log2FoldChange >= log2_cutoff ~ "Up",
    !is.na(x$padj) & x$padj <= padj_cutoff & x$log2FoldChange <= -log2_cutoff ~ "Down",
    TRUE ~ "NS"
  )

  out_file <- file.path(outdir, paste0("MA_", nm, ".pdf"))
  p <- ggplot2::ggplot(x, ggplot2::aes(log10(baseMean + 1), log2FoldChange, color = significance)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.65) +
    ggplot2::geom_hline(yintercept = c(-log2_cutoff, log2_cutoff), linetype = 2, color = "grey50") +
    ggplot2::scale_color_manual(values = c(Up = "#b2182b", Down = "#2166ac", NS = "grey75")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("MA plot:", nm),
      x = "log10(baseMean + 1)",
      y = "log2 fold change"
    )
  ggplot2::ggsave(out_file, plot = p, width = 5.5, height = 5)
  out_file
}

run_sig_heatmap <- function(vsd_mat, sig_genes, metadata, sample_column, group_column, outdir) {
  sig_mat <- vsd_mat[intersect(sig_genes, rownames(vsd_mat)), , drop = FALSE]
  if (nrow(sig_mat) < 2) return(NULL)
  if (!requireNamespace("pheatmap", quietly = TRUE)) return(NULL)
  ann_col <- metadata[, c(sample_column, group_column), drop = FALSE]
  rownames(ann_col) <- ann_col[[sample_column]]
  ann_col[[sample_column]] <- NULL
  out_file <- file.path(outdir, "significant_de_heatmap.pdf")
  pheatmap::pheatmap(
    scale_rows(sig_mat),
    show_rownames = FALSE,
    annotation_col = ann_col,
    filename = out_file,
    width = 7,
    height = 8
  )
  out_file
}
