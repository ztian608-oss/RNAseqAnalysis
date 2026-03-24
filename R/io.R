load_sample_table <- function(sample_table_file, count_dir) {
  assert_file_exists(sample_table_file, "sample_table")

  sample_table <- readr::read_csv(sample_table_file, show_col_types = FALSE) %>% as.data.frame()
  required_cols <- c("sampleName", "sampleFile", "group_name")
  if (!all(required_cols %in% colnames(sample_table))) {
    stop("sample_table must contain columns: sampleName, sampleFile, group_name")
  }

  if ("group_order" %in% colnames(sample_table)) {
    lv <- sample_table %>% dplyr::arrange(group_order) %>% dplyr::pull(group_name) %>% unique()
    sample_table$group_name <- factor(sample_table$group_name, levels = lv)
  } else {
    sample_table$group_name <- factor(sample_table$group_name)
  }

  count_paths <- file.path(count_dir, sample_table$sampleFile)
  missing_files <- sample_table$sampleFile[!file.exists(count_paths)]
  if (length(missing_files) > 0) {
    stop("Missing count files in count_dir: ", paste(missing_files, collapse = ", "))
  }

  sample_table
}

load_gene_info <- function(gtf_file) {
  assert_file_exists(gtf_file, "gtf")

  gtf <- rtracklayer::readGFF(gtf_file, version = 2L, tags = c("gene_id", "gene_name", "gene_type"))
  gene_info <- gtf %>%
    dplyr::select(gene_id, gene_name, gene_type) %>%
    unique() %>%
    dplyr::filter(!is.na(gene_id), !is.na(gene_name), gene_id != "", gene_name != "")

  if (nrow(gene_info) == 0) {
    stop("No valid gene annotation found in gtf (gene_id/gene_name).")
  }

  gene_info
}
