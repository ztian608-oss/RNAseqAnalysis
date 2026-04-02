assemble_result_table <- function(res, gene_info) {
  as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::left_join(gene_info, by = "gene_id") |>
    dplyr::relocate(gene_id, gene_name, gene_type)
}

extract_result_for_contrast <- function(dds, contrast_def, gene_info, variant) {
  cooks <- if (identical(variant$cooks_cutoff, "off")) FALSE else missing_arg()
  res <- if (identical(cooks, FALSE)) {
    DESeq2::results(
      dds,
      contrast = contrast_def,
      independentFiltering = variant$independent_filtering,
      alpha = variant$padj,
      cooksCutoff = FALSE
    )
  } else {
    DESeq2::results(
      dds,
      contrast = contrast_def,
      independentFiltering = variant$independent_filtering,
      alpha = variant$padj
    )
  }
  assemble_result_table(res, gene_info)
}

missing_arg <- function() {
  structure(list(), class = "missing_arg")
}
