normalize_organism_label <- function(organism) {
  x <- tolower(trimws(organism %||% "human"))
  if (x %in% c("human", "homo sapiens", "hs", "hsa")) return("human")
  if (x %in% c("mouse", "mus musculus", "mm", "mmu")) return("mouse")
  if (x %in% c("rat", "rattus norvegicus", "rn", "rno")) return("rat")
  x
}

canonicalize_symbols <- function(x) {
  tolower(trimws(as.character(x)))
}

interpretation_rules_path <- function() {
  root <- get0("skill_root", ifnotfound = NULL, inherits = TRUE)
  if (!is.null(root)) {
    return(file.path(root, "references", "interpretation_rules.json"))
  }
  file.path("references", "interpretation_rules.json")
}

load_interpretation_rules <- local({
  cache <- NULL
  function(force_reload = FALSE) {
    path <- interpretation_rules_path()
    if (!force_reload && !is.null(cache) && identical(cache$path, path)) {
      return(cache$data)
    }
    data <- jsonlite::read_json(path, simplifyVector = FALSE)
    cache <<- list(path = path, data = data)
    data
  }
})

infer_study_context_tags <- function(study_cfg) {
  text <- paste(
    study_cfg$objective %||% "",
    paste(study_cfg$disease_context %||% character(), collapse = " "),
    paste(study_cfg$focus_pathways %||% character(), collapse = " "),
    paste(study_cfg$focus_markers %||% character(), collapse = " ")
  )
  text <- tolower(text)

  tags <- character()
  rules <- load_interpretation_rules()
  for (rule in rules$context_tag_rules %||% list()) {
    pattern <- paste(rule$patterns %||% character(), collapse = "|")
    if (nzchar(pattern) && grepl(pattern, text)) {
      tags <- c(tags, rule$tag)
    }
  }
  unique(tags)
}

theme_library_generic <- function() {
  load_interpretation_rules()$generic_themes %||% list()
}

theme_library_context <- function() {
  load_interpretation_rules()$context_theme_sets %||% list()
}

theme_library_user_focus <- function(study_cfg) {
  markers <- unique(study_cfg$focus_markers %||% character())
  if (length(markers) == 0) return(list())
  list(
    user_priority_panel = list(
      label = "User Priority Panel",
      genes = markers
    )
  )
}

resolve_interpretation_rule_set <- function(study_cfg, organism = "human") {
  organism <- normalize_organism_label(organism)
  tags <- infer_study_context_tags(study_cfg)

  themes <- theme_library_generic()
  context_lib <- theme_library_context()
  for (tag in tags) {
    if (!is.null(context_lib[[tag]])) {
      themes <- c(themes, context_lib[[tag]])
    }
  }
  themes <- c(themes, theme_library_user_focus(study_cfg))

  list(
    organism = organism,
    context_tags = tags,
    themes = themes
  )
}

describe_interpretation_rule_set <- function(study_cfg, organism = "human") {
  rules <- resolve_interpretation_rule_set(study_cfg, organism)
  theme_labels <- vapply(rules$themes, function(x) x$label, character(1))
  c(
    paste0("- Interpretation organism context: ", rules$organism, "."),
    paste0("- Context tags inferred from the study description: ", if (length(rules$context_tags) == 0) "generic" else paste(rules$context_tags, collapse = ", "), "."),
    paste0("- Active interpretation panels: ", paste(theme_labels, collapse = ", "), ".")
  )
}

safe_gene_labels <- function(tbl) {
  labels <- if ("gene_name" %in% colnames(tbl)) tbl$gene_name else tbl$gene_id
  labels <- ifelse(is.na(labels) | labels == "", tbl$gene_id, labels)
  labels[!is.na(labels) & labels != ""]
}

sig_theme_hits <- function(tbl, variant_params, study_cfg, organism = "human", direction = c("up", "down")) {
  direction <- match.arg(direction)
  x <- tbl |>
    dplyr::filter(!is.na(.data$padj), !is.na(.data$log2FoldChange)) |>
    dplyr::filter(.data$padj <= variant_params$padj, abs(.data$log2FoldChange) >= variant_params$log2fc)

  if (direction == "up") {
    x <- x |> dplyr::filter(.data$log2FoldChange > 0)
  } else {
    x <- x |> dplyr::filter(.data$log2FoldChange < 0)
  }

  labels <- safe_gene_labels(x)
  label_key <- canonicalize_symbols(labels)
  rules <- resolve_interpretation_rule_set(study_cfg, organism)
  scores <- lapply(names(rules$themes), function(nm) {
    theme <- rules$themes[[nm]]
    theme_key <- canonicalize_symbols(theme$genes)
    hit_idx <- label_key %in% theme_key
    hits <- unique(labels[hit_idx])
    list(theme = nm, label = theme$label, n = length(hits), hits = hits)
  })
  names(scores) <- names(rules$themes)
  scores
}

top_theme_lines <- function(theme_scores, direction_label) {
  ranked <- theme_scores[order(vapply(theme_scores, `[[`, numeric(1), "n"), decreasing = TRUE)]
  ranked <- Filter(function(x) x$n > 0, ranked)
  if (length(ranked) == 0) {
    return(paste0("- No active interpretation panel was strongly represented among the ", direction_label, "-shift genes."))
  }

  top_ranked <- head(ranked, 2)
  vapply(top_ranked, function(x) {
    paste0(
      "- ", x$label,
      " signals are prominent in the ", direction_label, "-shift genes",
      " (examples: ", paste(head(x$hits, 5), collapse = ", "), ")."
    )
  }, character(1))
}

build_comparison_biological_summary <- function(tbl, comparison_name, study_cfg, variant_params, organism = "human") {
  up_scores <- sig_theme_hits(tbl, variant_params, study_cfg, organism, "up")
  down_scores <- sig_theme_hits(tbl, variant_params, study_cfg, organism, "down")

  lines <- c(
    paste0("- ", comparison_name, " should be interpreted at two levels: direction-specific DEG composition and behavior of the user-prioritized marker panel."),
    top_theme_lines(up_scores, "up"),
    top_theme_lines(down_scores, "down")
  )

  if (length(study_cfg$focus_pathways %||% character()) > 0) {
    lines <- c(lines, paste0("- These changes are most relevant to the requested pathways/processes: ", paste(study_cfg$focus_pathways, collapse = ", "), "."))
  }
  lines
}

cross_comparison_interpretation <- function(final_bundle, study_cfg) {
  marker_tbl <- final_bundle$focus_marker_table
  if (is.null(marker_tbl) || nrow(marker_tbl) == 0 || length(final_bundle$comparisons) < 2) {
    return("- Cross-comparison interpretation is limited because no shared focus-marker table is available.")
  }

  get_fc <- function(comparison, marker) {
    hit <- marker_tbl$log2FoldChange[marker_tbl$comparison == comparison &
      (marker_tbl$gene_name == marker | marker_tbl$gene_id == marker)]
    if (length(hit) == 0) return(NA_real_)
    hit[[1]]
  }

  comparisons <- final_bundle$comparisons
  comp_a <- comparisons[[1]]
  comp_b <- comparisons[[2]]
  markers <- unique(c("Tnnt2", "Actn2", "Myh6", "Myh7", "Atp2a2", "Ryr2", "Nppa", "Nppb", study_cfg$focus_markers))

  stronger_a <- character()
  stronger_b <- character()
  for (marker in unique(markers)) {
    a <- get_fc(comp_a, marker)
    b <- get_fc(comp_b, marker)
    if (!is.na(a) && !is.na(b) && abs(a - b) >= 0.5) {
      if (a > b) stronger_a <- c(stronger_a, marker)
      if (b > a) stronger_b <- c(stronger_b, marker)
    }
  }

  lines <- c("- Cross-comparison reading focuses on whether the two treatments push user-prioritized biological signals in the same direction and with similar magnitude.")
  if (length(stronger_a) > 0) {
    lines <- c(lines, paste0("- ", comp_a, " shows stronger relative shifts for markers such as ", paste(head(unique(stronger_a), 6), collapse = ", "), "."))
  }
  if (length(stronger_b) > 0) {
    lines <- c(lines, paste0("- ", comp_b, " shows stronger relative shifts for markers such as ", paste(head(unique(stronger_b), 6), collapse = ", "), "."))
  }
  if (length(stronger_a) == 0 && length(stronger_b) == 0) {
    lines <- c(lines, "- The treatment comparisons change the prioritized marker panel in largely similar magnitudes.")
  }
  lines
}
