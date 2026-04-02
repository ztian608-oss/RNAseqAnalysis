validate_numeric_counts <- function(count_mat) {
  if (!is.matrix(count_mat)) stop("count matrix must be a matrix")
  if (nrow(count_mat) == 0 || ncol(count_mat) == 0) stop("count matrix is empty")
  if (any(is.na(count_mat))) stop("count matrix contains NA")
  if (any(count_mat < 0)) stop("count matrix contains negative values")
}

validate_metadata_structure <- function(metadata, sample_column, group_column) {
  if (!(sample_column %in% colnames(metadata))) {
    stop("metadata is missing sample column: ", sample_column)
  }
  if (!(group_column %in% colnames(metadata))) {
    stop("metadata is missing group column: ", group_column)
  }
  if (anyDuplicated(metadata[[sample_column]]) > 0) {
    stop("metadata contains duplicated sample names in column: ", sample_column)
  }
  if (length(unique(as.character(metadata[[group_column]]))) < 2) {
    stop("At least two groups are required for differential analysis")
  }
}

align_counts_and_metadata <- function(count_mat, metadata, sample_column) {
  samples <- metadata[[sample_column]]
  if (!all(samples %in% colnames(count_mat))) {
    missing <- setdiff(samples, colnames(count_mat))
    stop("Samples missing from count matrix: ", paste(missing, collapse = ", "))
  }
  count_mat[, samples, drop = FALSE]
}

init_gene_info <- function(count_mat, gene_info) {
  if (is.null(gene_info)) {
    gene_info <- data.frame(
      gene_id = rownames(count_mat),
      gene_name = rownames(count_mat),
      gene_type = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  gene_info |> dplyr::filter(gene_id %in% rownames(count_mat))
}

prepare_analysis_input <- function(config) {
  metadata <- load_metadata_table(
    config$input$metadata,
    sample_column = config$design$sample_column,
    group_column = config$design$group_column
  )
  validate_metadata_structure(metadata, config$design$sample_column, config$design$group_column)

  if (config$input$counts_source == "matrix") {
    count_mat <- load_count_matrix(config$input$count_matrix)
  } else {
    sample_table <- load_sample_table_legacy(config$input$metadata, config$input$count_dir)
    config$design$sample_column <- "sampleName"
    config$design$group_column <- "group_name"
    metadata <- sample_table
    dds <- DESeq2::DESeqDataSetFromHTSeqCount(
      sampleTable = sample_table,
      directory = config$input$count_dir,
      design = as.formula(config$design$formula)
    )
    count_mat <- DESeq2::counts(dds)
  }

  validate_numeric_counts(count_mat)
  count_mat <- align_counts_and_metadata(count_mat, metadata, config$design$sample_column)
  gene_info <- init_gene_info(count_mat, load_annotation_table(config$input$annotation))
  background_genes <- load_background_genes(config$input$background_genes)

  list(
    config = config,
    metadata = metadata,
    count_mat = count_mat,
    gene_info = gene_info,
    background_genes = background_genes
  )
}
