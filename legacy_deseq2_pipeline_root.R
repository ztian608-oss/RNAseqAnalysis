#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(rtracklayer)
  library(tidyverse)
  library(pheatmap)
  library(EnhancedVolcano)
  library(SummarizedExperiment)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grep('^--file=', cmd_args)]
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub('^--file=', '', file_arg[[1]]))
} else {
  script_path <- normalizePath('legacy_deseq2_pipeline_root.R', mustWork = FALSE)
}
script_dir <- dirname(script_path)
module_files <- c(
  "legacy_cli_utils.R",
  "module_io.R",
  "legacy_analysis_helpers.R",
  "legacy_plot_helpers.R",
  "legacy_clustering_helpers.R"
)
for (mf in module_files) {
  full_path <- c(
    file.path(script_dir, "R", mf),
    file.path(script_dir, "scripts", "R", mf)
  )
  full_path <- full_path[file.exists(full_path)][1]
  if (!file.exists(full_path)) {
    stop("Module file not found for ", mf)
  }
  source(full_path)
}
timecourse_script <- c(
  file.path(script_dir, "R", "legacy_tcseq_timecourse_cli.R"),
  file.path(script_dir, "scripts", "R", "legacy_tcseq_timecourse_cli.R")
)
timecourse_script <- timecourse_script[file.exists(timecourse_script)][1]
if (!file.exists(timecourse_script)) {
  stop("Time-course clustering script not found.")
}
source(timecourse_script)

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  count_dir <- get_arg(args, "count_dir", required = TRUE)
  sample_table_file <- get_arg(args, "sample_table", required = TRUE)
  gtf_file <- get_arg(args, "gtf", required = TRUE)
  outdir <- get_arg(args, "outdir", default = "DESeq2_result")

  design_formula <- get_arg(args, "design", default = "~ group_name")
  ref_group <- get_arg(args, "ref_group", required = TRUE)
  contrast_file <- get_arg(args, "contrast_file", default = "")
  comparison_mode <- tolower(get_arg(args, "comparison_mode", default = "both"))

  padj_cutoff <- as.numeric(get_arg(args, "padj_cutoff", default = "0.01"))
  log2_cutoff <- as.numeric(get_arg(args, "log2_cutoff", default = "1"))
  min_count_mean <- as.numeric(get_arg(args, "min_count_mean", default = "10"))
  cluster_num <- as.integer(get_arg(args, "cluster_num", default = "6"))
  seed <- as.integer(get_arg(args, "seed", default = "12345"))
  run_volcano_plot <- to_logical(get_arg(args, "run_volcano", default = "TRUE"), TRUE)
  run_heatmap <- to_logical(get_arg(args, "run_heatmap", default = "TRUE"), TRUE)

  if (!comparison_mode %in% c("pairwise", "custom", "both")) {
    stop("comparison_mode must be one of: pairwise, custom, both")
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  sample_table <- load_sample_table_legacy(sample_table_file, count_dir)
  gene_info <- load_annotation_table(gtf_file)

  dds <- build_dds(
    sample_table = sample_table,
    count_dir = count_dir,
    design_formula = design_formula,
    gene_info = gene_info,
    min_count_mean = min_count_mean
  )

  vsd <- DESeq2::vst(dds, blind = FALSE)
  vsd_tbl <- as.data.frame(SummarizedExperiment::assay(vsd)) %>%
    tibble::rownames_to_column("gene_id") %>%
    tibble::as_tibble() %>%
    dplyr::right_join(gene_info, by = "gene_id") %>%
    dplyr::relocate(gene_id, gene_name, gene_type)
  readr::write_csv(vsd_tbl, file.path(outdir, "vsd.csv"))

  plot_deseq2_pca(
    vsd = vsd,
    output_pdf = file.path(outdir, "PCA_DESeq2.pdf"),
    output_csv = file.path(outdir, "PCA_DESeq2_data.csv"),
    intgroup = "group_name"
  )

  group_levels <- levels(SummarizedExperiment::colData(dds)$group_name)
  if (!(ref_group %in% group_levels)) {
    stop("ref_group is not in sample_table$group_name")
  }

  de_result_list <- list()
  if (comparison_mode %in% c("pairwise", "both")) {
    de_result_list <- c(de_result_list, run_pairwise_de(dds, group_levels, ref_group, gene_info, outdir))
  }

  if (comparison_mode %in% c("custom", "both")) {
    de_result_list <- c(de_result_list, run_custom_contrasts(dds, contrast_file, gene_info, outdir))
  }

  if (length(de_result_list) == 0) {
    stop("No DE comparison generated. Check comparison_mode/group labels/contrast_file.")
  }

  saveRDS(de_result_list, file.path(outdir, "de_result_list.rds"))

  sig_de_genes <- extract_significant_genes(de_result_list, log2_cutoff = log2_cutoff, padj_cutoff = padj_cutoff)
  readr::write_lines(sig_de_genes, file.path(outdir, "significant_de_genes.txt"))

  if (run_volcano_plot) {
    run_volcano(de_result_list, padj_cutoff, log2_cutoff, outdir)
  }

  vsd_mat <- SummarizedExperiment::assay(vsd)
  if (run_heatmap) {
    run_sig_heatmap(vsd_mat, sig_de_genes, sample_table, outdir)
  }

  if (length(group_levels) > 2) {
    run_timecourse_cluster(vsd_mat, sample_table, sig_de_genes, ref_group, cluster_num, outdir, seed = seed)
  } else {
    message("Group number <= 2; skip time-course clustering.")
  }

  message("Analysis finished. Results in: ", outdir)
}

main()
