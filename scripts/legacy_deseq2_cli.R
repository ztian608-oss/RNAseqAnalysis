#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
  cat(
    paste(
      "Legacy compatibility wrapper for RNAseqAnalysis.",
      "It converts the previous DESeq2-generic CLI arguments into the new config-driven workflow.",
      "",
      "Required:",
      "  --count_dir=PATH",
      "  --sample_table=FILE",
      "  --gtf=FILE",
      "  --ref_group=GROUP",
      "",
      "Optional:",
      "  --outdir=DIR",
      "  --design='~ group_name'",
      "  --contrast_file=FILE",
      "  --comparison_mode=pairwise|custom|both",
      "  --padj_cutoff=NUM",
      "  --log2_cutoff=NUM",
      "  --min_count_mean=NUM",
      "  --cluster_num=NUM",
      sep = "\n"
    )
  )
  quit(save = "no", status = 0)
}

parse_kv <- function(x) {
  kv <- strsplit(x, "=", fixed = TRUE)
  keys <- vapply(kv, function(y) sub("^--", "", y[[1]]), character(1))
  vals <- vapply(kv, function(y) if (length(y) > 1) paste(y[-1], collapse = "=") else "", character(1))
  as.list(stats::setNames(vals, keys))
}

arg <- parse_kv(args)
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])))
if (is.na(script_dir) || script_dir == "") {
  script_dir <- normalizePath("scripts", mustWork = FALSE)
}

module_files <- c(
  "module_logging.R", "module_config.R", "module_intake.R", "module_preflight.R",
  "module_io.R", "module_validation.R", "module_design.R",
  "module_deseq2_engine.R", "module_contrasts.R", "module_visualization.R", "module_deg_analysis.R",
  "module_tcseq_analysis.R", "module_enrichment.R", "module_stability.R", "module_result_selection.R",
  "module_report_generation.R", "module_workflow_orchestrator.R"
)
for (mf in module_files) {
  source(file.path(script_dir, "R", mf))
}

config_path <- file.path(tempdir(), "rnaseq_analysis_legacy_config.json")
cfg <- list(
  study = list(objective = "Legacy DESeq2 generic analysis"),
  input = list(
    counts_source = "htseq_dir",
    count_dir = arg$count_dir,
    metadata = arg$sample_table,
    annotation = arg$gtf,
    contrast_file = arg$contrast_file %||% ""
  ),
  design = list(
    formula = if (!is.null(arg$design) && nzchar(arg$design)) arg$design else "~ group_name",
    group_column = "group_name",
    sample_column = "sampleName",
    reference_group = arg$ref_group,
    comparison_mode = if (!is.null(arg$comparison_mode) && nzchar(arg$comparison_mode)) arg$comparison_mode else "pairwise"
  ),
  thresholds = list(
    padj = c(as.numeric(if (!is.null(arg$padj_cutoff) && nzchar(arg$padj_cutoff)) arg$padj_cutoff else "0.01")),
    log2fc = c(as.numeric(if (!is.null(arg$log2_cutoff) && nzchar(arg$log2_cutoff)) arg$log2_cutoff else "1")),
    min_count_mean = c(as.numeric(if (!is.null(arg$min_count_mean) && nzchar(arg$min_count_mean)) arg$min_count_mean else "10")),
    independent_filtering = c(TRUE),
    cooks_cutoff = c("default")
  ),
  trend = list(
    enabled = "auto",
    cluster_num = c(as.integer(if (!is.null(arg$cluster_num) && nzchar(arg$cluster_num)) arg$cluster_num else "6"))
  ),
  output = list(
    outdir = if (!is.null(arg$outdir) && nzchar(arg$outdir)) arg$outdir else "RNAseqAnalysis_output"
  )
)

jsonlite::write_json(cfg, config_path, auto_unbox = TRUE, pretty = TRUE)
config <- load_config(config_path)
show_preflight_plan(config)
result <- run_rnaseq_analysis_workflow(config)
cat(jsonlite::toJSON(result$manifest, auto_unbox = TRUE, pretty = TRUE), "\n")
