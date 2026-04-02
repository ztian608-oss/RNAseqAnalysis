aggregate_group_means <- function(vsd_mat, metadata, sample_column, group_column, ordered_groups) {
  group_means <- lapply(ordered_groups, function(g) {
    idx <- metadata[[sample_column]][metadata[[group_column]] == g]
    rowMeans(vsd_mat[, idx, drop = FALSE], na.rm = TRUE)
  })
  mat <- do.call(cbind, group_means)
  colnames(mat) <- ordered_groups
  mat
}

should_run_tcseq <- function(config, metadata) {
  group_count <- length(unique(as.character(metadata[[config$design$group_column]])))
  if (identical(config$trend$enabled, "true")) return(TRUE)
  if (identical(config$trend$enabled, "false")) return(FALSE)
  isTRUE(config$study$need_trend_analysis) || group_count >= config$trend$trigger_group_count_gt
}

run_tcseq_module <- function(vsd_mat, sig_genes, metadata, config, variant, variant_dir) {
  if (!should_run_tcseq(config, metadata) || length(sig_genes) < 2) {
    return(NULL)
  }
  if (!requireNamespace("TCseq", quietly = TRUE)) {
    return(list(
      out_dir = NULL,
      cluster_table = NULL,
      status = "TCseq package not installed"
    ))
  }
  groups <- levels(metadata[[config$design$group_column]])
  pattern <- aggregate_group_means(
    vsd_mat[sig_genes, , drop = FALSE],
    metadata,
    config$design$sample_column,
    config$design$group_column,
    groups
  )
  pattern <- scale_rows(pattern)
  pattern <- pattern[complete.cases(pattern), , drop = FALSE]
  if (nrow(pattern) < 2) return(NULL)

  tca <- TCseq::timeclust(pattern, algo = "cm", k = min(variant$cluster_num, nrow(pattern) - 1))
  cluster_tbl <- data.frame(
    gene_id = rownames(pattern),
    cluster = as.integer(tca@cluster),
    stringsAsFactors = FALSE
  )
  out_dir <- variant_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(cluster_tbl, file.path(out_dir, "cluster_assignment.csv"))

  grDevices::pdf(file.path(out_dir, "cluster_patterns.pdf"), width = 8, height = 6)
  p <- TCseq::timeclustplot(tca, cols = 2, title.size = 8, axis.text.size = 8, axis.title.size = 8)
  invisible(p)
  grDevices::dev.off()

  list(
    out_dir = out_dir,
    cluster_table = cluster_tbl,
    status = "completed"
  )
}
