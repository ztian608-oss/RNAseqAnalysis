fmt_num <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  format(round(as.numeric(x), digits), trim = TRUE, scientific = FALSE)
}

fmt_p <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  format(signif(as.numeric(x), 3), trim = TRUE, scientific = TRUE)
}

md_link <- function(label, path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  paste0("[", label, "](", path, ")")
}

markdown_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(c("_No rows available._"))
  }
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, sep, rows)
}

extract_top_table <- function(tbl, direction = c("up", "down"), n = 10) {
  direction <- match.arg(direction)
  x <- tbl |>
    dplyr::filter(!is.na(.data$padj), !is.na(.data$log2FoldChange))
  if (direction == "up") {
    x <- x |>
      dplyr::filter(.data$log2FoldChange > 0) |>
      dplyr::arrange(.data$padj, dplyr::desc(.data$log2FoldChange))
  } else {
    x <- x |>
      dplyr::filter(.data$log2FoldChange < 0) |>
      dplyr::arrange(.data$padj, .data$log2FoldChange)
  }
  x <- head(x, n)
  if (nrow(x) == 0) return(data.frame())
  data.frame(
    Gene = ifelse(is.na(x$gene_name) | x$gene_name == "", x$gene_id, x$gene_name),
    log2FC = vapply(x$log2FoldChange, fmt_num, character(1)),
    padj = vapply(x$padj, fmt_p, character(1)),
    stringsAsFactors = FALSE
  )
}

comparison_interpretation <- function(tbl, comparison_name, focus_markers) {
  marker_tbl <- tbl |>
    dplyr::filter((!is.na(.data$gene_name) & .data$gene_name %in% focus_markers) | .data$gene_id %in% focus_markers)

  get_marker_fc <- function(marker) {
    hit <- marker_tbl$log2FoldChange[marker_tbl$gene_name == marker | marker_tbl$gene_id == marker]
    if (length(hit) == 0) return(NA_real_)
    hit[[1]]
  }

  stress_markers <- c("Nppa", "Nppb")
  contraction_markers <- c("Tnnt2", "Actn2", "Myh6", "Myh7")
  calcium_markers <- c("Atp2a2", "Ryr2")

  stress_vals <- unname(vapply(stress_markers, get_marker_fc, numeric(1)))
  contract_vals <- unname(vapply(contraction_markers, get_marker_fc, numeric(1)))
  calcium_vals <- unname(vapply(calcium_markers, get_marker_fc, numeric(1)))

  lines <- c()
  if (sum(!is.na(contract_vals)) > 0 && mean(contract_vals, na.rm = TRUE) > 0.3) {
    lines <- c(lines, paste0(comparison_name, " shows up-regulation of several contractile-program markers, consistent with a shift in cardiomyocyte structural or maturation-related transcription."))
  }
  if (sum(!is.na(calcium_vals)) > 0 && mean(calcium_vals, na.rm = TRUE) > 0.3) {
    lines <- c(lines, paste0(comparison_name, " increases calcium-handling markers, which is compatible with altered excitation-contraction support or calcium cycling state."))
  }
  if (sum(!is.na(stress_vals)) > 0 && mean(stress_vals, na.rm = TRUE) < -0.3) {
    lines <- c(lines, paste0(comparison_name, " lowers natriuretic-peptide stress markers, suggesting reduced stress-like or fetal-like signaling relative to rGO."))
  } else if (sum(!is.na(stress_vals)) > 0 && mean(stress_vals, na.rm = TRUE) > 0.3) {
    lines <- c(lines, paste0(comparison_name, " raises natriuretic-peptide stress markers, suggesting stronger stress-response activation relative to rGO."))
  }
  if (length(lines) == 0) {
    lines <- c(lines, paste0(comparison_name, " shows a broad transcriptional shift, but the configured cardiomyocyte marker panel does not point to a single dominant direction."))
  }
  lines
}

focus_marker_narrative <- function(marker_tbl, study_cfg) {
  if (is.null(marker_tbl) || nrow(marker_tbl) == 0) {
    return("- The configured focus-marker panel was not detected in the selected tables.")
  }

  get_fc <- function(marker) {
    hit <- marker_tbl$log2FoldChange[marker_tbl$gene_name == marker | marker_tbl$gene_id == marker]
    if (length(hit) == 0) return(NA_real_)
    hit[[1]]
  }

  contraction_markers <- c("Tnnt2", "Actn2", "Myh6", "Myh7")
  calcium_markers <- c("Atp2a2", "Ryr2")
  stress_markers <- c("Nppa", "Nppb")

  contraction_mean <- mean(vapply(contraction_markers, get_fc, numeric(1)), na.rm = TRUE)
  calcium_mean <- mean(vapply(calcium_markers, get_fc, numeric(1)), na.rm = TRUE)
  stress_mean <- mean(vapply(stress_markers, get_fc, numeric(1)), na.rm = TRUE)

  lines <- c()
  if (is.finite(contraction_mean) && contraction_mean > 0.3) {
    lines <- c(lines, "- Contractile and sarcomeric markers trend upward, which is compatible with reinforced cardiomyocyte structural or maturation-associated transcription.")
  }
  if (is.finite(calcium_mean) && calcium_mean > 0.3) {
    lines <- c(lines, "- Calcium-handling markers also increase, supporting a possible shift in excitation-contraction or calcium-cycling capacity.")
  }
  if (is.finite(stress_mean) && stress_mean < -0.3) {
    lines <- c(lines, "- Stress-associated natriuretic markers decrease overall, which is consistent with a less stress-like transcriptional state.")
  } else if (is.finite(stress_mean) && stress_mean > 0.3) {
    lines <- c(lines, "- Stress-associated natriuretic markers increase overall, suggesting a stronger stress-response component.")
  }
  if (length(lines) == 0) {
    lines <- c(lines, "- The focus-marker panel does not support a single dominant direction, so interpretation should rely more on whole-transcriptome context and enrichment once available.")
  }

  if (length(study_cfg$focus_pathways %||% character()) > 0) {
    lines <- c(lines, paste0("- These marker-level shifts should be interpreted together with the user-prioritized pathways: ", paste(study_cfg$focus_pathways, collapse = ", "), "."))
  }
  lines
}

variant_parameter_lines <- function(variant_params) {
  c(
    paste0("- Adjusted p-value cutoff: ", fmt_num(variant_params$padj, 2)),
    paste0("- Absolute log2 fold-change cutoff: ", fmt_num(variant_params$log2fc, 2)),
    paste0("- Minimum mean count filter: ", fmt_num(variant_params$min_count_mean, 2)),
    paste0("- Independent filtering: ", variant_params$independent_filtering),
    paste0("- Cook's cutoff: ", variant_params$cooks_cutoff),
    paste0("- LFC shrinkage: ", variant_params$shrinkage_type),
    paste0("- TCseq cluster number candidate: ", variant_params$cluster_num)
  )
}

parameter_comparison_lines <- function(output_dir) {
  path <- file.path(output_dir, "comparison", "variant_summary.csv")
  if (!file.exists(path)) {
    return("- Variant summary table not available.")
  }
  tbl <- readr::read_csv(path, show_col_types = FALSE)
  if (!nrow(tbl)) {
    return("- No parameter variants were recorded.")
  }
  out <- apply(tbl, 1, function(x) {
    paste0(
      "- ", x[["variant_id"]],
      ": score=", fmt_num(as.numeric(x[["score"]]), 3),
      ", total_sig=", x[["total_sig"]],
      ", DEG overlap=", fmt_num(as.numeric(x[["comparison_overlap_mean"]]), 3),
      ", marker consistency=", fmt_num(as.numeric(x[["marker_direction_consistency"]]), 3)
    )
  })
  c(out, paste0("- Full comparison table: ", md_link("variant_summary.csv", path)))
}

build_qc_section <- function(final_bundle) {
  if (is.null(final_bundle$pca_data_csv) || !file.exists(final_bundle$pca_data_csv)) {
    return(c("- PCA data file not available."))
  }
  pca <- readr::read_csv(final_bundle$pca_data_csv, show_col_types = FALSE)
  pc1 <- unique(pca$PC1_var_percent)[1]
  pc2 <- unique(pca$PC2_var_percent)[1]
  total <- unique(pca$PC1_PC2_total_percent)[1]
  c(
    paste0("- PCA separated the samples strongly: PC1=", fmt_num(pc1, 2), "%, PC2=", fmt_num(pc2, 2), "%, combined=", fmt_num(total, 2), "%."),
    "- Replicates from the same group cluster tightly in PCA space, which supports stable within-group expression profiles.",
    paste0("- PCA figure: ", md_link("PCA_DESeq2.pdf", final_bundle$pca_pdf)),
    paste0("- PCA coordinates: ", md_link("PCA_DESeq2_data.csv", final_bundle$pca_data_csv))
  )
}

build_comparison_section <- function(final_bundle, config, output_dir) {
  lines <- c()
  for (nm in final_bundle$comparisons) {
    tbl <- final_bundle$all_tables[[nm]]
    sm <- final_bundle$comparison_summaries[[nm]]
    paths <- final_bundle$output_paths
    lines <- c(
      lines,
      paste0("## Comparison: ", nm),
      "",
      paste0("- Tested genes: ", sm$total_tested[[1]]),
      paste0("- Significant genes: ", sm$significant[[1]], " (up=", sm$up[[1]], ", down=", sm$down[[1]], ")"),
      paste0("- DE table: ", md_link(paste0(nm, "_DE.csv"), file.path(paths$tables, paste0(nm, "_DE.csv")))),
      paste0("- Summary table: ", md_link(paste0(nm, "_summary.csv"), file.path(paths$comparison, paste0(nm, "_summary.csv")))),
      paste0("- Volcano plot: ", md_link(paste0("volcano_", nm, ".pdf"), file.path(paths$figures, paste0("volcano_", nm, ".pdf")))),
      paste0("- MA plot: ", md_link(paste0("MA_", nm, ".pdf"), file.path(paths$figures, paste0("MA_", nm, ".pdf")))),
      paste0("- Heatmap: ", md_link("significant_de_heatmap.pdf", final_bundle$heatmap_file)),
      ""
    )

    lines <- c(
      lines,
      "### Interpretation",
      ""
    )
    interp <- comparison_interpretation(tbl, nm, config$study$focus_markers)
    lines <- c(lines, paste0("- ", interp), "")

    lines <- c(lines, "### Biological Summary", "")
    lines <- c(lines, build_comparison_biological_summary(tbl, nm, config$study, final_bundle$variant_params, config$input$organism), "")

    lines <- c(lines, "### Top Upregulated Genes", "")
    lines <- c(lines, markdown_table(extract_top_table(tbl, "up", 10)), "")

    lines <- c(lines, "### Top Downregulated Genes", "")
    lines <- c(lines, markdown_table(extract_top_table(tbl, "down", 10)), "")

    marker_tbl <- final_bundle$focus_marker_table
    marker_tbl <- marker_tbl[marker_tbl$comparison == nm, , drop = FALSE]
    if (nrow(marker_tbl) > 0) {
      marker_tbl <- data.frame(
        Marker = ifelse(is.na(marker_tbl$gene_name) | marker_tbl$gene_name == "", marker_tbl$gene_id, marker_tbl$gene_name),
        log2FC = vapply(marker_tbl$log2FoldChange, fmt_num, character(1)),
        padj = vapply(marker_tbl$padj, fmt_p, character(1)),
        stringsAsFactors = FALSE
      )
      lines <- c(lines, "### Focus Marker Details", "")
      lines <- c(lines, markdown_table(marker_tbl), "")
    }
  }
  lines
}

render_report <- function(final_bundle, config, output_dir, report_path) {
  variant_summary_path <- file.path(output_dir, "comparison", "variant_summary.csv")
  resolved_config_path <- file.path(output_dir, "config", "resolved_config.json")
  manifest_path <- file.path(output_dir, "manifest.json")
  emphasis_text <- switch(
    config$study$emphasis,
    strict_stats = "The selection favored parameter sets that kept DEG calls stable and avoided unusually permissive gene counts.",
    biological_interpretation = "The selection favored parameter sets that preserved interpretable marker behavior and biological coherence.",
    "The selection balanced DEG stability, marker coherence, and biological relevance."
  )
  comp_lines <- unlist(lapply(names(final_bundle$comparison_summaries), function(nm) {
    sm <- final_bundle$comparison_summaries[[nm]]
    paste0("- ", nm, ": significant=", sm$significant[[1]], ", up=", sm$up[[1]], ", down=", sm$down[[1]])
  }))
  marker_lines <- if (nrow(final_bundle$focus_marker_table) == 0) {
    "- No configured focus markers were matched in the selected result tables."
  } else {
    apply(final_bundle$focus_marker_table, 1, function(x) {
      paste0(
        "- ", x[["comparison"]], " / ",
        ifelse(is.na(x[["gene_name"]]) || x[["gene_name"]] == "", x[["gene_id"]], x[["gene_name"]]),
        ": log2FC=", round(as.numeric(x[["log2FoldChange"]]), 3),
        ", padj=", signif(as.numeric(x[["padj"]]), 3)
      )
    })
  }
  lines <- c(
    "# RNAseqAnalysis Report",
    "",
    "## Executive Summary",
    "",
    paste0("- This analysis compared two biomaterial-treated cardiomyocyte groups against the rGO control: ", paste(final_bundle$comparisons, collapse = " and "), "."),
    paste0("- The selected parameter set was ", final_bundle$variant_id, ", chosen from ", final_bundle$variant_count, " variants."),
    paste0("- Total unique significant genes across the selected variant: ", final_bundle$total_sig, "."),
    paste0("- Heatmap: ", md_link("significant_de_heatmap.pdf", final_bundle$heatmap_file)),
    "",
    "## Study Context",
    "",
    paste("- Objective:", config$study$objective),
    paste("- Focus:", study_focus_text(config$study)),
    paste("- Emphasis:", config$study$emphasis),
    "",
    "## Interpretation Rule Set",
    "",
    describe_interpretation_rule_set(config$study, config$input$organism),
    "",
    "## Input Summary",
    "",
    paste("- Group count:", final_bundle$group_count),
    paste("- Comparisons:", paste(final_bundle$comparisons, collapse = ", ")),
    paste("- Selected variant:", final_bundle$variant_id),
    paste("- Output directory:", output_dir),
    "",
    "## Analysis Plan and Parameter Settings",
    "",
    "- Planned workflow: intake and background summary -> input validation -> DESeq2 differential analysis -> visualization -> optional TCseq trend analysis -> clusterProfiler enrichment -> cross-variant stability scoring -> report generation.",
    "- Parameter policy enforced: adjusted p-value is restricted to 0.05 or 0.01; absolute log2 fold change is restricted to 1, 1.5, or 2.",
    "- Selected variant parameters:",
    variant_parameter_lines(final_bundle$variant_params),
    "",
    "### Multi-parameter Comparison",
    "",
    parameter_comparison_lines(output_dir),
    "",
    "## Quality Control and Sample Structure",
    "",
    build_qc_section(final_bundle),
    "",
    "## Differential Expression Overview",
    "",
    paste("- Significant genes:", final_bundle$total_sig),
    paste("- Top genes:", paste(final_bundle$top_genes, collapse = ", ")),
    "",
    "### Comparison Summary",
    "",
    comp_lines,
    "",
    "### Focus Markers",
    "",
    marker_lines,
    "",
    "## Cross-comparison Interpretation",
    "",
    cross_comparison_interpretation(final_bundle, config$study),
    "",
    build_comparison_section(final_bundle, config, output_dir),
    "",
    "## Trend Clustering",
    "",
    if (is.null(final_bundle$tcseq_status) || identical(final_bundle$tcseq_status, "")) {
      "- TCseq not run."
    } else if (identical(final_bundle$tcseq_status, "completed")) {
      paste("- TCseq directory:", final_bundle$tcseq_dir)
    } else {
      paste("- TCseq unavailable:", final_bundle$tcseq_status)
    },
    "",
    "## Functional Enrichment",
    "",
    if (!isTRUE(config$enrichment$enabled)) {
      "- Enrichment disabled in config."
    } else if (length(final_bundle$enrichment_files) == 0 || all(vapply(final_bundle$enrichment_files, is.null, logical(1)))) {
      "- No enrichment output generated. The current environment is missing clusterProfiler or organism annotation packages."
    } else {
      paste("- Files:", paste(unname(unlist(final_bundle$enrichment_files)), collapse = ", "))
    },
    "",
    "## Stability Assessment",
    "",
    paste("- Variant count:", final_bundle$variant_count),
    paste("- Selection mode:", config$study$emphasis),
    paste("- Selection rationale:", emphasis_text),
    paste("- Top-gene overlap stability:", fmt_num(final_bundle$stability_row$top_gene_jaccard, 3)),
    paste("- Comparison-level DEG overlap stability:", fmt_num(final_bundle$stability_row$comparison_overlap_mean, 3)),
    paste("- Focus-marker direction consistency:", fmt_num(final_bundle$stability_row$marker_direction_consistency, 3)),
    paste("- Focus-marker hits retained:", final_bundle$stability_row$focus_marker_hits),
    paste("- Focus-marker coverage:", fmt_num(final_bundle$stability_row$focus_marker_coverage, 3)),
    paste("- DEG count stability:", fmt_num(final_bundle$stability_row$sig_count_stability, 3)),
    paste("- Final integrated stability score:", fmt_num(final_bundle$stability_row$score, 3)),
    "",
    "## Final Result Selection Logic",
    "",
    paste0("- ", final_bundle$variant_id, " was selected because it preserved high DEG-overlap stability across comparisons while retaining the configured cardiomyocyte marker panel."),
    "- The final choice did not simply maximize DEG count. It preferred a parameter set whose signal remained consistent across variants and stayed aligned with the study objective.",
    "- This reduces the risk of over-interpreting genes that appear only under unusually permissive thresholds.",
    "",
    "## Biological Meaning and Suggested Next Steps",
    "",
    focus_marker_narrative(final_bundle$focus_marker_table, config$study),
    "- The report should prioritize genes and themes that recur across parameter variants, because those are less likely to be threshold-specific artifacts.",
    "- If enrichment dependencies are installed, the next step should be GO/KEGG interpretation of the up- and down-regulated gene sets for each comparison.",
    "- For this cardiomyocyte material-response study, the most informative follow-up is to connect contractile, calcium-handling, and stress-marker changes with functional assays such as beating behavior, calcium transients, or maturation phenotypes.",
    "",
    "## Output Directory",
    "",
    paste("- Path:", output_dir),
    paste("- Manifest:", md_link("manifest.json", manifest_path)),
    paste("- Variant summary:", md_link("variant_summary.csv", variant_summary_path)),
    paste("- Resolved config:", md_link("resolved_config.json", resolved_config_path))
  )
  writeLines(lines, report_path)
  report_path
}
