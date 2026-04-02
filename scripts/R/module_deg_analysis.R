apply_shrinkage <- function(dds, tbl, contrast_name, variant) {
  if (identical(variant$shrinkage_type, "none")) {
    return(tbl)
  }
  rn <- DESeq2::resultsNames(dds)
  coef_name <- rn[grepl(contrast_name, rn, fixed = TRUE)]
  if (length(coef_name) == 0) {
    return(tbl)
  }
  shr <- tryCatch(
    DESeq2::lfcShrink(dds, coef = coef_name[[1]], type = variant$shrinkage_type, quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(shr)) return(tbl)
  shr_tbl <- assemble_result_table(shr, tbl[, c("gene_id", "gene_name", "gene_type")])
  cols <- intersect(c("gene_id", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), colnames(shr_tbl))
  dplyr::left_join(
    dplyr::select(tbl, -dplyr::any_of(c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))),
    dplyr::select(shr_tbl, dplyr::all_of(cols)),
    by = "gene_id"
  )
}

summarize_deg_table <- function(tbl, variant) {
  sig <- tbl |>
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
    dplyr::filter(padj <= variant$padj, abs(log2FoldChange) >= variant$log2fc)

  list(
    significant = sig,
    summary = data.frame(
      total_tested = sum(!is.na(tbl$padj)),
      significant = nrow(sig),
      up = sum(sig$log2FoldChange > 0),
      down = sum(sig$log2FoldChange < 0)
    )
  )
}

ranked_gene_list <- function(tbl, id_col = "gene_id") {
  x <- tbl |>
    dplyr::filter(!is.na(stat)) |>
    dplyr::arrange(dplyr::desc(stat))
  stats::setNames(x$stat, x[[id_col]])
}

top_feature_labels <- function(tbl, n = 20) {
  x <- tbl |>
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange)))
  labels <- if ("gene_name" %in% colnames(x)) x$gene_name else x$gene_id
  labels <- ifelse(is.na(labels) | labels == "", x$gene_id, labels)
  labels <- labels[!is.na(labels) & labels != ""]
  unique(head(labels, n))
}

extract_focus_marker_table <- function(all_tables, focus_markers) {
  if (length(focus_markers) == 0) {
    return(data.frame())
  }
  rows <- list()
  for (nm in names(all_tables)) {
    tbl <- all_tables[[nm]]
    hit <- tbl |>
      dplyr::filter((!is.na(.data$gene_name) & .data$gene_name %in% focus_markers) | .data$gene_id %in% focus_markers) |>
      dplyr::mutate(comparison = nm) |>
      dplyr::select("comparison", "gene_id", "gene_name", "log2FoldChange", "padj")
    rows[[nm]] <- hit
  }
  dplyr::bind_rows(rows)
}
