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
  for (nm in names(de_result_list)) {
    tbl <- as.data.frame(de_result_list[[nm]])
    tbl <- tbl[!is.na(tbl$padj) & !is.na(tbl$log2FoldChange), ]
    if (nrow(tbl) == 0) next

    up_top <- tbl %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::slice_min(padj, n = 5) %>% dplyr::pull(gene_name)
    down_top <- tbl %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::slice_min(padj, n = 5) %>% dplyr::pull(gene_name)

    p <- EnhancedVolcano::EnhancedVolcano(
      tbl,
      lab = tbl$gene_name,
      selectLab = tbl$gene_name[tbl$gene_name %in% c(up_top, down_top)],
      x = "log2FoldChange",
      y = "padj",
      subtitle = nm,
      pCutoff = padj_cutoff,
      FCcutoff = log2_cutoff,
      drawConnectors = TRUE,
      widthConnectors = 0.2,
      colConnectors = "gray50",
      labSize = 3
    )

    ggplot2::ggsave(file.path(outdir, paste0("volcanoplot_", nm, ".pdf")), plot = p, width = 5, height = 5)
  }
}
