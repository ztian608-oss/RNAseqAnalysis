normalize_intake <- function(study_cfg) {
  study_cfg$objective <- study_cfg$objective %||% "Run differential expression analysis."
  study_cfg$disease_context <- unname(unlist(study_cfg$disease_context %||% character()))
  study_cfg$focus_pathways <- unname(unlist(study_cfg$focus_pathways %||% character()))
  study_cfg$focus_markers <- unname(unlist(study_cfg$focus_markers %||% character()))
  study_cfg$focus_groups <- unname(unlist(study_cfg$focus_groups %||% character()))
  study_cfg$emphasis <- study_cfg$emphasis %||% "balanced"
  study_cfg$need_trend_analysis <- isTRUE(study_cfg$need_trend_analysis)
  study_cfg
}

study_focus_text <- function(study_cfg) {
  parts <- c(
    study_cfg$objective,
    if (length(study_cfg$disease_context) > 0) paste("Disease:", paste(study_cfg$disease_context, collapse = ", ")) else NULL,
    if (length(study_cfg$focus_pathways) > 0) paste("Pathways:", paste(study_cfg$focus_pathways, collapse = ", ")) else NULL,
    if (length(study_cfg$focus_markers) > 0) paste("Markers:", paste(study_cfg$focus_markers, collapse = ", ")) else NULL
  )
  paste(parts, collapse = " | ")
}
