plot_deseq2_pca <- function(vsd, output_pdf, output_csv, intgroup = "group_name") {
  pca_data <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"), 2)
  pc12_total <- round(sum(percent_var[1:2]), 2)

  pca_data$PC1_var_percent <- percent_var[1]
  pca_data$PC2_var_percent <- percent_var[2]
  pca_data$PC1_PC2_total_percent <- pc12_total

  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(PC1, PC2, color = .data[[intgroup]], label = name)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "PCA (DESeq2 VST)",
      subtitle = paste0("PC1+PC2 explained variance: ", pc12_total, "%"),
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = intgroup,
      caption = "Variance values are derived from DESeq2::plotPCA percentVar"
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave(output_pdf, plot = p, width = 5, height = 4)
  readr::write_csv(pca_data, output_csv)
}

run_volcano <- function(de_result_list, padj_cutoff, log2_cutoff, outdir) {
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    stop("EnhancedVolcano package is required for volcano plot generation.")
  }
  for (nm in names(de_result_list)) {
    tbl <- as.data.frame(de_result_list[[nm]])
    tbl <- tbl[!is.na(tbl$pvalue) & !is.na(tbl$log2FoldChange), ]
    if (nrow(tbl) == 0) next

    labels <- ifelse(is.na(tbl$gene_name) | tbl$gene_name == "", tbl$gene_id, tbl$gene_name)
    up_top <- tbl %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::slice_min(padj, n = 5) %>% dplyr::pull(gene_name)
    down_top <- tbl %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::slice_min(padj, n = 5) %>% dplyr::pull(gene_name)
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
      selectLab = labels[labels %in% c(up_top, down_top)],
      x = "log2FoldChange",
      y = "pvalue",
      title = paste("Volcano plot:", nm),
      subtitle = paste0("x = log2FoldChange; y = -log10(p value); point classes use padj <= ", padj_cutoff,
        " and |log2FC| >= ", log2_cutoff),
      xlab = bquote(~Log[2]~ "fold change"),
      ylab = bquote(~-Log[10]~ italic(P)),
      pCutoff = 1,
      FCcutoff = 0,
      drawConnectors = TRUE,
      widthConnectors = 0.2,
      colConnectors = "gray50",
      labSize = 3,
      pointSize = 1.8,
      colAlpha = 0.9,
      colCustom = keyvals
    )

    ggplot2::ggsave(file.path(outdir, paste0("volcanoplot_", nm, ".pdf")), plot = p, width = 5, height = 5)
  }
}
