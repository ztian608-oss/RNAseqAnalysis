has_study_background <- function(study_cfg) {
  any(nzchar(c(
    study_cfg$objective %||% "",
    paste(study_cfg$disease_context %||% character(), collapse = ""),
    paste(study_cfg$focus_pathways %||% character(), collapse = ""),
    paste(study_cfg$focus_markers %||% character(), collapse = "")
  )))
}

build_analysis_plan <- function(config) {
  steps <- c(
    "1. Validate count input, metadata, annotation, and contrast settings.",
    "2. Normalize design formula and reference group.",
    "3. Expand parameter variants and run DESeq2 differential analysis.",
    "4. Generate PCA, DEG tables, volcano plots, and DEG heatmap.",
    if (identical(config$trend$enabled, "false")) NULL else "5. Run TCseq trend clustering when group count and config require it.",
    "6. Run clusterProfiler enrichment on DEG sets and TCseq clusters when available.",
    "7. Compare parameter variants and select the final result set.",
    "8. Write manifest, variant summary, and markdown report."
  )
  unique(steps)
}

missing_background_questions <- function() {
  c(
    "What disease, biological system, or treatment context does this study focus on?",
    "Are there pathways or marker genes that matter most for interpretation?",
    "Should the final result favor stricter statistics or broader biological interpretability?"
  )
}

show_preflight_plan <- function(config) {
  cat("Preflight plan for RNAseqAnalysis\n")
  cat(paste(build_analysis_plan(config), collapse = "\n"), "\n")
  if (!has_study_background(config$study)) {
    cat("\nMissing background context. Recommended questions before trusting final variant selection:\n")
    cat(paste0("- ", missing_background_questions(), collapse = "\n"), "\n")
  }
}
