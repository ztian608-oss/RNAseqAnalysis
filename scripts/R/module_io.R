assert_file_exists <- function(path, desc) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(desc, " not found: ", path)
  }
}

load_metadata_table <- function(metadata_file, sample_column, group_column) {
  assert_file_exists(metadata_file, "metadata")
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE) |> as.data.frame()
  required_cols <- c(sample_column, group_column)
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  metadata[[group_column]] <- factor(metadata[[group_column]])
  metadata
}

load_count_matrix <- function(matrix_file, sample_column = "gene_id") {
  assert_file_exists(matrix_file, "count_matrix")
  tbl <- readr::read_csv(matrix_file, show_col_types = FALSE)
  gene_col <- if (sample_column %in% colnames(tbl)) sample_column else colnames(tbl)[1]
  mat <- as.data.frame(tbl)
  rownames(mat) <- mat[[gene_col]]
  mat[[gene_col]] <- NULL
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

load_sample_table_legacy <- function(sample_table_file, count_dir) {
  assert_file_exists(sample_table_file, "sample_table")
  sample_table <- readr::read_csv(sample_table_file, show_col_types = FALSE) |> as.data.frame()
  required_cols <- c("sampleName", "sampleFile", "group_name")
  missing_cols <- setdiff(required_cols, colnames(sample_table))
  if (length(missing_cols) > 0) {
    stop("sample_table must contain columns: sampleName, sampleFile, group_name")
  }
  sample_table$group_name <- factor(sample_table$group_name)
  count_paths <- file.path(count_dir, sample_table$sampleFile)
  missing_files <- sample_table$sampleFile[!file.exists(count_paths)]
  if (length(missing_files) > 0) {
    stop("Missing count files in count_dir: ", paste(missing_files, collapse = ", "))
  }
  sample_table
}

load_annotation_table <- function(annotation_file) {
  if (is.null(annotation_file) || !nzchar(annotation_file)) {
    return(NULL)
  }
  assert_file_exists(annotation_file, "annotation")
  ext <- tools::file_ext(annotation_file)
  if (tolower(ext) %in% c("gtf", "gff", "gff3")) {
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gtf <- rtracklayer::readGFF(annotation_file, version = 2L, tags = c("gene_id", "gene_name", "gene_type"))
      gene_info <- gtf |>
        dplyr::select(gene_id, gene_name, gene_type) |>
        unique()
    } else {
      raw_gtf <- readr::read_tsv(
        annotation_file,
        comment = "#",
        col_names = FALSE,
        show_col_types = FALSE,
        progress = FALSE
      )
      if (ncol(raw_gtf) < 9) {
        stop("GTF/GFF must contain at least 9 columns")
      }
      attr_col <- raw_gtf[[9]]
      extract_attr <- function(x, key) {
        out <- stringr::str_match(x, paste0(key, ' "?([^\";]+)"?'))[, 2]
        out[is.na(out)] <- ""
        out
      }
      gene_info <- tibble::tibble(
        gene_id = extract_attr(attr_col, "gene_id"),
        gene_name = extract_attr(attr_col, "gene_name"),
        gene_type = extract_attr(attr_col, "gene_type")
      ) |>
        dplyr::filter(.data$gene_id != "") |>
        unique()
    }
  } else {
    gene_info <- readr::read_csv(annotation_file, show_col_types = FALSE) |>
      dplyr::rename_with(~tolower(.x))
    if (!("gene_id" %in% colnames(gene_info))) {
      gene_info <- dplyr::rename(gene_info, gene_id = 1)
    }
    if (!("gene_name" %in% colnames(gene_info))) {
      gene_info$gene_name <- gene_info$gene_id
    }
    if (!("gene_type" %in% colnames(gene_info))) {
      gene_info$gene_type <- NA_character_
    }
    gene_info <- dplyr::select(gene_info, gene_id, gene_name, gene_type)
  }
  gene_info |>
    dplyr::filter(!is.na(gene_id), gene_id != "") |>
    unique()
}

load_background_genes <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  assert_file_exists(path, "background_genes")
  readr::read_lines(path)
}
