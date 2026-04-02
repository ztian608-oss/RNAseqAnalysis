`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

default_config <- function() {
  list(
    study = list(
      objective = "Run differential expression analysis.",
      disease_context = character(),
      focus_pathways = character(),
      focus_markers = character(),
      focus_groups = character(),
      emphasis = "balanced",
      need_trend_analysis = FALSE
    ),
    input = list(
      counts_source = "matrix",
      metadata = NULL,
      annotation = NULL,
      contrast_file = "",
      background_genes = "",
      organism = "human",
      id_type = "SYMBOL"
    ),
    design = list(
      formula = "~ group_name",
      group_column = "group_name",
      sample_column = "sampleName",
      reference_group = NULL,
      comparison_mode = "pairwise",
      run_lrt = FALSE,
      lrt_reduced = ""
    ),
    thresholds = list(
      padj = c(0.05),
      log2fc = c(1),
      min_count_mean = c(10),
      independent_filtering = c(TRUE),
      cooks_cutoff = c("default")
    ),
    shrinkage = list(
      enabled = TRUE,
      type = c("apeglm")
    ),
    trend = list(
      enabled = "auto",
      trigger_group_count_gt = 3L,
      cluster_num = c(4L, 6L),
      seed = 12345L
    ),
    enrichment = list(
      enabled = TRUE,
      mode = c("ORA"),
      databases = c("GO_BP", "KEGG"),
      universe = "expressed_genes",
      custom_term2gene = ""
    ),
    output = list(
      outdir = "RNAseqAnalysis_output",
      report_name = "RNAseqAnalysis_report.md"
    )
  )
}

merge_config <- function(base, override) {
  for (nm in names(override)) {
    if (is.list(base[[nm]]) && is.list(override[[nm]])) {
      base[[nm]] <- merge_config(base[[nm]], override[[nm]])
    } else {
      base[[nm]] <- override[[nm]]
    }
  }
  base
}

load_config <- function(path) {
  cfg <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  normalize_config(merge_config(default_config(), cfg), config_path = path)
}

normalize_config <- function(config, config_path = NULL) {
  config$meta <- list(
    config_path = config_path %||% "",
    loaded_at = timestamp_now()
  )
  config$study <- normalize_intake(config$study)
  config <- normalize_parameter_policy(config)
  config
}

normalize_parameter_policy <- function(config) {
  config$thresholds$padj <- unique(as.numeric(unlist(config$thresholds$padj %||% c(0.05))))
  config$thresholds$log2fc <- unique(as.numeric(unlist(config$thresholds$log2fc %||% c(1))))
  config$thresholds$min_count_mean <- unique(as.numeric(unlist(config$thresholds$min_count_mean %||% c(10))))
  config$thresholds$independent_filtering <- unique(as.logical(unlist(config$thresholds$independent_filtering %||% c(TRUE))))
  config$thresholds$cooks_cutoff <- unique(as.character(unlist(config$thresholds$cooks_cutoff %||% c("default"))))
  config$trend$cluster_num <- unique(as.integer(unlist(config$trend$cluster_num %||% c(4L, 6L))))
  validate_parameter_policy(config)
}

validate_parameter_policy <- function(config) {
  allowed_padj <- c(0.05, 0.01)
  allowed_log2fc <- c(1, 1.5, 2)
  allowed_clusters <- c(4L, 6L, 8L)

  bad_padj <- setdiff(config$thresholds$padj, allowed_padj)
  if (length(bad_padj) > 0) {
    stop("Invalid padj thresholds: ", paste(bad_padj, collapse = ", "),
         ". Allowed values are 0.05 and 0.01.")
  }

  bad_lfc <- setdiff(config$thresholds$log2fc, allowed_log2fc)
  if (length(bad_lfc) > 0) {
    stop("Invalid log2FC thresholds: ", paste(bad_lfc, collapse = ", "),
         ". Allowed values are 1, 1.5, and 2.")
  }

  bad_clusters <- setdiff(config$trend$cluster_num, allowed_clusters)
  if (length(bad_clusters) > 0) {
    stop("Invalid TCseq cluster numbers: ", paste(bad_clusters, collapse = ", "),
         ". Allowed values are 4, 6, and 8.")
  }

  bad_min_count <- config$thresholds$min_count_mean[!is.finite(config$thresholds$min_count_mean) | config$thresholds$min_count_mean <= 0]
  if (length(bad_min_count) > 0) {
    stop("Invalid min_count_mean thresholds: ", paste(bad_min_count, collapse = ", "),
         ". Values must be > 0, and 10 is the recommended default.")
  }

  config
}

expand_variants <- function(config) {
  grid <- expand.grid(
    padj = config$thresholds$padj,
    log2fc = config$thresholds$log2fc,
    min_count_mean = config$thresholds$min_count_mean,
    independent_filtering = config$thresholds$independent_filtering,
    cooks_cutoff = config$thresholds$cooks_cutoff,
    shrinkage_type = if (isTRUE(config$shrinkage$enabled)) config$shrinkage$type else "none",
    cluster_num = config$trend$cluster_num,
    stringsAsFactors = FALSE
  )

  variants <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    variants[[i]] <- list(
      variant_id = sprintf("variant_%02d", i),
      padj = as.numeric(row$padj[[1]]),
      log2fc = as.numeric(row$log2fc[[1]]),
      min_count_mean = as.numeric(row$min_count_mean[[1]]),
      independent_filtering = as.logical(row$independent_filtering[[1]]),
      cooks_cutoff = row$cooks_cutoff[[1]],
      shrinkage_type = row$shrinkage_type[[1]],
      cluster_num = as.integer(row$cluster_num[[1]])
    )
  }
  variants
}
