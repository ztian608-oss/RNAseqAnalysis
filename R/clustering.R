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

run_timecourse_cluster <- function(vsd_mat, sample_table, sig_genes, cluster_num, outdir, seed = 12345) {
  if (length(sig_genes) == 0) {
    message("No significant DE genes found; skip time-course clustering.")
    return(invisible(NULL))
  }

  group_levels <- levels(sample_table$group_name)
  means_by_group <- lapply(group_levels, function(g) {
    idx <- sample_table$sampleName[sample_table$group_name == g]
    rowMeans(vsd_mat[sig_genes, idx, drop = FALSE], na.rm = TRUE)
  })

  pattern_mat <- do.call(cbind, means_by_group)
  colnames(pattern_mat) <- group_levels
  pattern_scaled <- scale_rows(pattern_mat)

  set.seed(seed)
  km <- kmeans(pattern_scaled, centers = cluster_num, nstart = 50)

  ann_row <- data.frame(cluster = factor(km$cluster))
  rownames(ann_row) <- rownames(pattern_scaled)

  pheatmap::pheatmap(
    pattern_scaled[order(km$cluster), , drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    annotation_row = ann_row[order(km$cluster), , drop = FALSE],
    filename = file.path(outdir, "timecourse_cluster_heatmap.pdf"),
    width = 6,
    height = 8
  )

  cluster_df <- tibble::tibble(
    gene_id = rownames(pattern_scaled),
    cluster = km$cluster
  )
  readr::write_csv(cluster_df, file.path(outdir, "timecourse_cluster_result.csv"))

  centers_long <- as.data.frame(km$centers) %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(-cluster, names_to = "group_name", values_to = "z_expr")

  p <- ggplot2::ggplot(
    centers_long,
    ggplot2::aes(x = group_name, y = z_expr, group = cluster, color = cluster)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Time-course cluster patterns", x = "Group", y = "Scaled expression")

  ggplot2::ggsave(file.path(outdir, "timecourse_cluster_patterns.pdf"), plot = p, width = 7, height = 4)
}
