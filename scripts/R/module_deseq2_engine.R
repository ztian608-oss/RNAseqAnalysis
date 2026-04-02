build_dds_from_matrix <- function(count_mat, metadata, design_formula, min_count_mean) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = S4Vectors::DataFrame(metadata),
    design = as.formula(design_formula)
  )
  dds <- dds[MatrixGenerics::rowMeans(DESeq2::counts(dds)) >= min_count_mean, ]
  if (nrow(dds) == 0) {
    stop("No genes remain after min_count_mean filtering")
  }
  dds
}

run_deseq_fit <- function(dds, config) {
  if (isTRUE(config$design$run_lrt)) {
    reduced <- config$design$lrt_reduced
    if (is.null(reduced) || !nzchar(reduced)) {
      stop("LRT requires design$lrt_reduced")
    }
    DESeq2::DESeq(dds, test = "LRT", reduced = as.formula(reduced), quiet = TRUE)
  } else {
    DESeq2::DESeq(dds, quiet = TRUE)
  }
}
