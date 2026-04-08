write_manifest <- function(manifest, outdir) {
  manifest_path <- file.path(outdir, "manifest.json")
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  manifest_path
}

initialize_output_layout <- function(outdir) {
  dirs <- list(
    root = outdir,
    config = file.path(outdir, "config"),
    deseq2 = file.path(outdir, "deseq2"),
    tcseq = file.path(outdir, "tcseq"),
    enrichment = file.path(outdir, "enrichment"),
    figures = file.path(outdir, "figures"),
    tables = file.path(outdir, "tables"),
    comparison = file.path(outdir, "comparison"),
    report = file.path(outdir, "report"),
    logs = file.path(outdir, "logs")
  )
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}

# Shared output paths for non-TCseq analyses (only one copy regardless of variants)
shared_output_paths <- function(layout) {
  list(
    deseq2 = layout$deseq2,
    enrichment = layout$enrichment,
    figures = layout$figures,
    tables = layout$tables,
    comparison = layout$comparison
  )
}

# TCseq-specific paths (per variant since clustering results differ)
tcseq_output_paths <- function(layout, variant_id) {
  list(
    tcseq = file.path(layout$tcseq, variant_id),
    figures = file.path(layout$figures, paste0("tcseq_", variant_id))
  )
}

# Legacy function for compatibility
variant_output_paths <- function(layout, variant_id) {
  list(
    deseq2 = file.path(layout$deseq2, variant_id),
    tcseq = file.path(layout$tcseq, variant_id),
    enrichment = file.path(layout$enrichment, variant_id),
    figures = file.path(layout$figures, variant_id),
    tables = file.path(layout$tables, variant_id),
    comparison = file.path(layout$comparison, variant_id)
  )
}

write_run_log <- function(layout, lines) {
  log_path <- file.path(layout$logs, "run.log")
  writeLines(lines, log_path)
  log_path
}

# Check if variants differ only in TCseq parameters (cluster_num)
variants_differ_only_in_tcseq <- function(variants) {
  if (length(variants) <= 1) return(FALSE)
  # Check if all variants have same DE-related parameters
  for (i in 2:length(variants)) {
    v1 <- variants[[1]]
    vi <- variants[[i]]
    if (v1$padj != vi$padj || v1$log2fc != vi$log2fc || 
        v1$min_count_mean != vi$min_count_mean || 
        v1$independent_filtering != vi$independent_filtering ||
        v1$cooks_cutoff != vi$cooks_cutoff) {
      return(FALSE)
    }
  }
  TRUE
}

run_rnaseq_analysis_workflow <- function(config) {
  layout <- initialize_output_layout(config$output$outdir)
  input <- prepare_analysis_input(config)
  config$design$reference_group <- resolve_reference_group(
    input$metadata,
    config$design$group_column,
    config$design$reference_group
  )
  input$metadata <- normalize_group_factor(input$metadata, config$design$group_column, config$design$reference_group)
  variants <- expand_variants(config)
  comparisons <- load_contrasts(config, input$metadata)
  group_count <- length(unique(as.character(input$metadata[[config$design$group_column]])))
  run_log <- c(
    paste("Run started:", timestamp_now()),
    paste("Output root:", config$output$outdir),
    paste("Comparisons:", paste(names(comparisons), collapse = ", ")),
    paste("Group count:", group_count)
  )

  # Determine if we can share non-TCseq outputs across variants
  share_outputs <- variants_differ_only_in_tcseq(variants)
  
  if (share_outputs) {
    message("[Info] Variants differ only in TCseq cluster_num. ",
            "DE analysis and enrichment will run once, TCseq will run per variant.")
  }

  # Create shared paths for non-TCseq outputs
  shared <- shared_output_paths(layout)
  for (d in c(shared$deseq2, shared$enrichment, shared$figures, shared$tables, shared$comparison)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  variant_results <- vector("list", length(variants))
  
  # Use first variant for DE analysis if sharing
  base_variant <- variants[[1]]
  
  # ============================================
  # PHASE 1: Run DE analysis, enrichment, figures (once if sharing)
  # ============================================
  message("[Progress] Running differential expression analysis...")
  run_log <- c(run_log, paste("DE analysis parameters:",
    "padj=", base_variant$padj,
    "log2fc=", base_variant$log2fc,
    "min_count_mean=", base_variant$min_count_mean))
  
  dds_bundle <- build_dds_from_matrix(
    input$count_mat,
    input$metadata,
    config$design$formula,
    base_variant$min_count_mean
  )
  dds <- dds_bundle$dds
  dds <- run_deseq_fit(dds, config)
  vsd <- DESeq2::vst(dds, blind = FALSE)
  vsd_tbl <- as.data.frame(SummarizedExperiment::assay(vsd)) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::left_join(input$gene_info, by = "gene_id")
  readr::write_csv(vsd_tbl, file.path(shared$deseq2, "vsd.csv"))
  
  # PCA plot
  plot_deseq2_pca(
    vsd,
    file.path(shared$figures, "PCA_DESeq2.pdf"),
    file.path(shared$tables, "PCA_DESeq2_data.csv"),
    config$design$group_column
  )

  # Comparison loop
  comp_summaries <- list()
  all_sig <- character()
  all_tables <- list()
  significant_gene_sets <- list()
  enrich_files <- list()
  enrichment_plot_files <- list()
  
  for (nm in names(comparisons)) {
    message("[Progress] Processing comparison: ", nm)
    tbl <- extract_result_for_contrast(dds, comparisons[[nm]], input$gene_info, base_variant)
    tbl <- apply_shrinkage(dds, tbl, nm, base_variant)
    deg <- summarize_deg_table(tbl, base_variant)
    readr::write_csv(tbl, file.path(shared$tables, paste0(nm, "_DE.csv")))
    readr::write_csv(deg$summary, file.path(shared$comparison, paste0(nm, "_summary.csv")))
    run_volcano(tbl, nm, base_variant$padj, base_variant$log2fc, shared$figures)
    run_ma_plot(tbl, nm, base_variant$padj, base_variant$log2fc, shared$figures)
    comp_summaries[[nm]] <- deg$summary
    all_tables[[nm]] <- tbl
    significant_gene_sets[[nm]] <- unique(deg$significant$gene_id)
    all_sig <- unique(c(all_sig, deg$significant$gene_id))
    up_genes <- unique(deg$significant$gene_id[deg$significant$log2FoldChange > 0])
    down_genes <- unique(deg$significant$gene_id[deg$significant$log2FoldChange < 0])
    enrich_files[[nm]] <- list(
      up = run_ora_for_gene_set(
        up_genes,
        config,
        resolve_universe(config, input$background_genes, input$gene_info, tbl),
        paste0(nm, "_up"),
        shared$enrichment
      ),
      down = run_ora_for_gene_set(
        down_genes,
        config,
        resolve_universe(config, input$background_genes, input$gene_info, tbl),
        paste0(nm, "_down"),
        shared$enrichment
      )
    )
    enrichment_plot_files[[nm]] <- run_enrichment_bubble_plots(
      enrich_files[[nm]],
      comparison_name = comparison_label_from_name(nm, comparisons[[nm]]),
      outdir = shared$figures
    )
  }

  # Heatmap
  heatmap_file <- run_sig_heatmap(
    SummarizedExperiment::assay(vsd),
    all_sig,
    input$metadata,
    config$design$sample_column,
    config$design$group_column,
    shared$figures
  )

  # ============================================
  # PHASE 2: Run TCseq per variant (only if needed)
  # ============================================
  group_count_for_trend <- length(unique(as.character(input$metadata[[config$design$group_column]])))
  need_trend <- isTRUE(config$trend$enabled) || 
                (identical(config$trend$enabled, "auto") && group_count_for_trend > config$trend$trigger_group_count_gt)
  
  if (need_trend && length(variants) > 0) {
    message("[Progress] Running TCseq trend analysis for ", length(variants), " variant(s)...")
    for (i in seq_along(variants)) {
      variant <- variants[[i]]
      tcseq_paths <- tcseq_output_paths(layout, variant$variant_id)
      dir.create(tcseq_paths$tcseq, recursive = TRUE, showWarnings = FALSE)
      dir.create(tcseq_paths$figures, recursive = TRUE, showWarnings = FALSE)
      
      message("[Progress] TCseq variant ", i, "/", length(variants), 
              ": cluster_num=", variant$cluster_num)
      run_log <- c(run_log, paste("TCseq:", variant$variant_id, 
        "cluster_num=", variant$cluster_num))
      
      tcseq_res <- run_tcseq_module(
        SummarizedExperiment::assay(vsd),
        all_sig,
        input$metadata,
        config,
        variant,
        tcseq_paths$tcseq
      )
      
      # Store variant result with TCseq info
      variant_results[[i]] <- list(
        variant_id = variant$variant_id,
        variant_params = variant,
        filter_summary = dds_bundle$filter_summary,
        total_sig = length(all_sig),
        top_genes = head(unique(unlist(lapply(all_tables, top_feature_labels, n = 10))), 20),
        comparisons = names(comparisons),
        contrast_defs = comparisons,
        all_tables = all_tables,
        significant_gene_sets = significant_gene_sets,
        comparison_summaries = comp_summaries,
        focus_marker_table = extract_focus_marker_table(all_tables, config$study$focus_markers),
        output_paths = list(
          deseq2 = shared$deseq2,
          enrichment = shared$enrichment,
          figures = shared$figures,
          tables = shared$tables,
          comparison = shared$comparison,
          tcseq = tcseq_paths$tcseq
        ),
        vsd_csv = file.path(shared$deseq2, "vsd.csv"),
        pca_pdf = file.path(shared$figures, "PCA_DESeq2.pdf"),
        pca_data_csv = file.path(shared$tables, "PCA_DESeq2_data.csv"),
        heatmap_file = heatmap_file,
        tcseq_dir = if (is.null(tcseq_res)) NULL else tcseq_res$out_dir,
        tcseq_status = if (is.null(tcseq_res)) "" else tcseq_res$status %||% "",
        tcseq_figures_dir = tcseq_paths$figures,
        enrichment_files = enrich_files,
        enrichment_plot_files = enrichment_plot_files
      )
    }
  } else {
    # No TCseq needed - just record variant info without TCseq dirs
    for (i in seq_along(variants)) {
      variant <- variants[[i]]
      variant_results[[i]] <- list(
        variant_id = variant$variant_id,
        variant_params = variant,
        filter_summary = dds_bundle$filter_summary,
        total_sig = length(all_sig),
        top_genes = head(unique(unlist(lapply(all_tables, top_feature_labels, n = 10))), 20),
        comparisons = names(comparisons),
        contrast_defs = comparisons,
        all_tables = all_tables,
        significant_gene_sets = significant_gene_sets,
        comparison_summaries = comp_summaries,
        focus_marker_table = extract_focus_marker_table(all_tables, config$study$focus_markers),
        output_paths = shared,
        vsd_csv = file.path(shared$deseq2, "vsd.csv"),
        pca_pdf = file.path(shared$figures, "PCA_DESeq2.pdf"),
        pca_data_csv = file.path(shared$tables, "PCA_DESeq2_data.csv"),
        heatmap_file = heatmap_file,
        tcseq_dir = NULL,
        tcseq_status = "not_run",
        tcseq_figures_dir = NULL,
        enrichment_files = enrich_files,
        enrichment_plot_files = enrichment_plot_files
      )
    }
  }

  stability_tbl <- compare_variants(variant_results, config$study)
  variant_summary_path <- file.path(layout$comparison, "variant_summary.csv")
  readr::write_csv(stability_tbl, variant_summary_path)
  selected <- select_final_variant(stability_tbl)
  final_bundle <- variant_results[[match(selected$variant_id, vapply(variant_results, `[[`, character(1), "variant_id"))]]
  final_bundle$variant_count <- length(variant_results)
  final_bundle$group_count <- group_count
  final_bundle$stability_row <- selected

  config_path_out <- file.path(layout$config, "resolved_config.json")
  jsonlite::write_json(config, config_path_out, auto_unbox = TRUE, pretty = TRUE)
  log_path <- write_run_log(layout, run_log)
  report_path <- file.path(layout$report, config$output$report_name)
  render_report(final_bundle, config, config$output$outdir, report_path)

  manifest <- list(
    outdir = config$output$outdir,
    layout = layout,
    final_variant = final_bundle$variant_id,
    group_count = final_bundle$group_count,
    comparisons = final_bundle$comparisons,
    stability_summary = list(
      variant_count = nrow(stability_tbl),
      selected_reason = paste("Selected", final_bundle$variant_id, "based on DEG overlap stability, marker consistency, and study emphasis")
    ),
    files = list(
      report_md = report_path,
      config_json = config_path_out,
      manifest_json = file.path(config$output$outdir, "manifest.json"),
      variant_summary_csv = variant_summary_path,
      pca_pdf = final_bundle$pca_pdf,
      pca_data_csv = final_bundle$pca_data_csv,
      vsd_csv = final_bundle$vsd_csv,
      selected_deseq2_dir = final_bundle$output_paths$deseq2,
      selected_figures_dir = final_bundle$output_paths$figures,
      selected_tables_dir = final_bundle$output_paths$tables,
      selected_comparison_dir = final_bundle$output_paths$comparison,
      selected_enrichment_dir = final_bundle$output_paths$enrichment,
      selected_tcseq_dir = final_bundle$output_paths$tcseq,
      log_txt = log_path
    )
  )
  manifest_path <- write_manifest(manifest, config$output$outdir)
  manifest$files$manifest_json <- manifest_path

  list(manifest = manifest, final_bundle = final_bundle, stability = stability_tbl)
}
