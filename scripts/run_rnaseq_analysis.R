#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])))
if (is.na(script_dir) || script_dir == "") {
  script_dir <- normalizePath("scripts", mustWork = FALSE)
}
skill_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

module_files <- c(
  "module_logging.R", "module_config.R", "module_intake.R", "module_preflight.R",
  "module_io.R", "module_validation.R", "module_design.R",
  "module_deseq2_engine.R", "module_contrasts.R", "module_visualization.R", "module_deg_analysis.R",
  "module_biological_interpretation.R",
  "module_tcseq_analysis.R", "module_enrichment.R", "module_stability.R", "module_result_selection.R",
  "module_report_generation.R", "module_workflow_orchestrator.R"
)

for (mf in module_files) {
  source(file.path(script_dir, "R", mf))
}

args <- commandArgs(trailingOnly = TRUE)
config_arg <- args[grepl("^--config=", args)]
if (length(config_arg) == 0) {
  stop("Usage: Rscript scripts/run_rnaseq_analysis.R --config=path/to/config.json")
}
config_path <- sub("^--config=", "", config_arg[[1]])

config <- load_config(config_path)
show_preflight_plan(config)
result <- run_rnaseq_analysis_workflow(config)
cat(jsonlite::toJSON(result$manifest, auto_unbox = TRUE, pretty = TRUE), "\n")
