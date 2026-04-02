resolve_reference_group <- function(metadata, group_column, configured_reference = NULL) {
  groups <- unique(as.character(metadata[[group_column]]))
  if (!is.null(configured_reference) && nzchar(configured_reference)) {
    if (!(configured_reference %in% groups)) {
      stop("reference_group not found in metadata: ", configured_reference)
    }
    return(configured_reference)
  }
  groups[[1]]
}

normalize_group_factor <- function(metadata, group_column, ref_group) {
  metadata[[group_column]] <- factor(as.character(metadata[[group_column]]))
  metadata[[group_column]] <- stats::relevel(metadata[[group_column]], ref = ref_group)
  metadata
}

load_contrasts <- function(config, metadata) {
  mode <- config$design$comparison_mode %||% "pairwise"
  groups <- levels(metadata[[config$design$group_column]])
  ref_group <- config$design$reference_group
  out <- list()

  if (mode %in% c("pairwise", "both")) {
    for (g in groups[groups != ref_group]) {
      out[[paste0(g, "_vs_", ref_group)]] <- c(config$design$group_column, g, ref_group)
    }
  }

  if (mode %in% c("custom", "both")) {
    assert_file_exists(config$input$contrast_file, "contrast_file")
    tbl <- readr::read_csv(config$input$contrast_file, show_col_types = FALSE)
    for (i in seq_len(nrow(tbl))) {
      nm <- paste0(tbl$treat[[i]], "_vs_", tbl$base[[i]])
      out[[nm]] <- c(config$design$group_column, tbl$treat[[i]], tbl$base[[i]])
    }
  }

  if (length(out) == 0) {
    stop("No contrasts were generated")
  }
  out
}
