run_sig_heatmap <- function(vsd_mat, sig_genes, sample_table, outdir) {
  if (length(sig_genes) == 0) {
    message("No significant DE genes found; skip significant DE heatmap.")
    return(invisible(NULL))
  }

  sig_mat <- vsd_mat[intersect(sig_genes, rownames(vsd_mat)), , drop = FALSE]
  if (nrow(sig_mat) == 0) {
    message("No overlap between significant genes and VST matrix; skip heatmap.")
    return(invisible(NULL))
  }

  sig_scaled <- scale_rows(sig_mat)
  ann_col <- data.frame(group_name = sample_table$group_name)
  rownames(ann_col) <- sample_table$sampleName

  pheatmap::pheatmap(
    sig_scaled,
    show_rownames = FALSE,
    clustering_method = "ward.D2",
    annotation_col = ann_col,
    filename = file.path(outdir, "significant_de_heatmap.pdf"),
    width = 7,
    height = 8
  )
}

# Legacy TCseq helper moved to scripts/R/legacy_tcseq_timecourse_cli.R
