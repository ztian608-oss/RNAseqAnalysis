#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(rtracklayer)
  library(tidyverse)
  library(pheatmap)
  library(EnhancedVolcano)
})

parse_args <- function(args) {
  if (length(args) == 0) {
    stop("No arguments provided. Use --key=value format.")
  }
  kv <- strsplit(args, "=", fixed = TRUE)
  keys <- vapply(kv, function(x) gsub("^--", "", x[[1]]), character(1))
  vals <- vapply(kv, function(x) {
    if (length(x) < 2) "" else paste(x[-1], collapse = "=")
  }, character(1))
  as.list(setNames(vals, keys))
}

get_arg <- function(arg_list, key, default = NULL, required = FALSE) {
  val <- arg_list[[key]]
  if (is.null(val) || val == "") {
    if (required) {
      stop(paste0("Missing required argument: --", key))
    }
    return(default)
  }
  val
}

to_logical <- function(x, default = FALSE) {
  if (is.null(x) || x == "") return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  s[s == 0] <- 1
  (x - m) / s
}

plot_deseq2_pca <- function(vsd, output_pdf, output_csv, intgroup = "group_name") {
  pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(PC1, PC2, color = .data[[intgroup]], label = name)) +
    geom_point(size = 3) +
    labs(
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance"),
      color = intgroup
    ) +
    theme_bw()

  ggsave(output_pdf, plot = p, width = 5, height = 4)
  write_csv(pca_data, output_csv)
}

run_pairwise_de <- function(dds, group_levels, ref_group, gene_info, outdir) {
  de_result_list <- list()

  cmp_list <- group_levels[group_levels != ref_group]
  for (cmp in cmp_list) {
    dds_sub <- dds[, dds$group_name %in% c(ref_group, cmp)]
    dds_sub$group_name <- droplevels(dds_sub$group_name)
    dds_sub <- DESeq(dds_sub)

    res <- results(dds_sub, contrast = c("group_name", cmp, ref_group))
    vsd_sub <- vst(dds_sub, blind = FALSE)

    res_tbl <- cbind(as.data.frame(assay(vsd_sub)), as.data.frame(res)) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      right_join(gene_info, by = "gene_id") %>%
      relocate(gene_id, gene_name, gene_type)

    file <- file.path(outdir, paste0(cmp, "_vs_", ref_group, "_DE.csv"))
    write_csv(res_tbl, file)
    de_result_list[[paste0(cmp, "_vs_", ref_group)]] <- res_tbl
  }

  de_result_list
}

run_custom_contrasts <- function(dds, contrast_file, gene_info, outdir) {
  if (is.null(contrast_file) || contrast_file == "") return(list())

  contrast_tbl <- read_csv(contrast_file, show_col_types = FALSE)
  required_cols <- c("treat", "base")
  if (!all(required_cols %in% colnames(contrast_tbl))) {
    stop("contrast_file must contain columns: treat, base")
  }

  out <- list()
  for (i in seq_len(nrow(contrast_tbl))) {
    treat <- contrast_tbl$treat[[i]]
    base <- contrast_tbl$base[[i]]

    dds_sub <- dds[, dds$group_name %in% c(base, treat)]
    dds_sub$group_name <- droplevels(dds_sub$group_name)
    dds_sub <- DESeq(dds_sub)

    res <- results(dds_sub, contrast = c("group_name", treat, base))
    vsd_sub <- vst(dds_sub, blind = FALSE)

    res_tbl <- cbind(as.data.frame(assay(vsd_sub)), as.data.frame(res)) %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      right_join(gene_info, by = "gene_id") %>%
      relocate(gene_id, gene_name, gene_type)

    nm <- paste0(treat, "_vs_", base)
    write_csv(res_tbl, file.path(outdir, paste0(nm, "_DE.csv")))
    out[[nm]] <- res_tbl
  }

  out
}

run_volcano <- function(de_result_list, padj_cutoff, log2_cutoff, outdir) {
  for (nm in names(de_result_list)) {
    tbl <- as.data.frame(de_result_list[[nm]])
    tbl <- tbl[!is.na(tbl$padj) & !is.na(tbl$log2FoldChange), ]
    if (nrow(tbl) == 0) next

    up_top <- tbl %>% filter(log2FoldChange > 0) %>% slice_min(padj, n = 5) %>% pull(gene_name)
    down_top <- tbl %>% filter(log2FoldChange < 0) %>% slice_min(padj, n = 5) %>% pull(gene_name)

    p <- EnhancedVolcano(
      tbl,
      lab = tbl$gene_name,
      selectLab = tbl$gene_name[tbl$gene_name %in% c(up_top, down_top)],
      x = "log2FoldChange",
      y = "padj",
      subtitle = nm,
      pCutoff = padj_cutoff,
      FCcutoff = log2_cutoff,
      drawConnectors = TRUE,
      widthConnectors = 0.2,
      colConnectors = "gray50",
      labSize = 3
    )

    ggsave(file.path(outdir, paste0("volcanoplot_", nm, ".pdf")), plot = p, width = 5, height = 5)
  }
}

run_timecourse_cluster <- function(vsd_mat, sample_table, sig_genes, cluster_num, outdir) {
  if (length(sig_genes) == 0) {
    message("No significant DE genes found; skip time-course clustering.")
    return(invisible(NULL))
  }

  group_levels <- levels(sample_table$group_name)
  means_by_group <- lapply(group_levels, function(g) {
    idx <- sample_table$sampleName[sample_table$group_name == g]
    rowMeans(vsd_mat[sig_genes, idx, drop = FALSE], na.rm = TRUE)
  })

  pattern_mat <- do.call(cbind, means_by_group)
  colnames(pattern_mat) <- group_levels
  pattern_scaled <- scale_rows(pattern_mat)

  set.seed(12345)
  km <- kmeans(pattern_scaled, centers = cluster_num, nstart = 50)

  ann_row <- data.frame(cluster = factor(km$cluster))
  rownames(ann_row) <- rownames(pattern_scaled)

  pheatmap(
    pattern_scaled[order(km$cluster), , drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    annotation_row = ann_row[order(km$cluster), , drop = FALSE],
    filename = file.path(outdir, "timecourse_cluster_heatmap.pdf"),
    width = 6,
    height = 8
  )

  cluster_df <- tibble(
    gene_id = rownames(pattern_scaled),
    cluster = km$cluster
  )
  write_csv(cluster_df, file.path(outdir, "timecourse_cluster_result.csv"))

  centers_long <- as.data.frame(km$centers) %>%
    rownames_to_column("cluster") %>%
    pivot_longer(-cluster, names_to = "group_name", values_to = "z_expr")

  p <- ggplot(centers_long, aes(x = group_name, y = z_expr, group = cluster, color = cluster)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_bw() +
    labs(title = "Time-course cluster patterns", x = "Group", y = "Scaled expression")

  ggsave(file.path(outdir, "timecourse_cluster_patterns.pdf"), plot = p, width = 7, height = 4)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  count_dir <- get_arg(args, "count_dir", required = TRUE)
  sample_table_file <- get_arg(args, "sample_table", required = TRUE)
  gtf_file <- get_arg(args, "gtf", required = TRUE)
  outdir <- get_arg(args, "outdir", default = "DESeq2_result")

  design_formula <- get_arg(args, "design", default = "~ group_name")
  ref_group <- get_arg(args, "ref_group", required = TRUE)
  contrast_file <- get_arg(args, "contrast_file", default = "")

  padj_cutoff <- as.numeric(get_arg(args, "padj_cutoff", default = "0.01"))
  log2_cutoff <- as.numeric(get_arg(args, "log2_cutoff", default = "1"))
  min_count_mean <- as.numeric(get_arg(args, "min_count_mean", default = "10"))
  cluster_num <- as.integer(get_arg(args, "cluster_num", default = "6"))
  seed <- as.integer(get_arg(args, "seed", default = "12345"))
  run_volcano_plot <- to_logical(get_arg(args, "run_volcano", default = "TRUE"), TRUE)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  sample_table <- read_csv(sample_table_file, show_col_types = FALSE) %>% as.data.frame()
  required_cols <- c("sampleName", "sampleFile", "group_name")
  if (!all(required_cols %in% colnames(sample_table))) {
    stop("sample_table must contain columns: sampleName, sampleFile, group_name")
  }

  if ("group_order" %in% colnames(sample_table)) {
    lv <- sample_table %>% arrange(group_order) %>% pull(group_name) %>% unique()
    sample_table$group_name <- factor(sample_table$group_name, levels = lv)
  } else {
    sample_table$group_name <- factor(sample_table$group_name)
  }

  gtf <- readGFF(gtf_file, version = 2L, tags = c("gene_id", "gene_name", "gene_type"))
  gene_info <- gtf %>%
    dplyr::select(gene_id, gene_name, gene_type) %>%
    unique() %>%
    filter(!is.na(gene_name))

  design <- as.formula(design_formula)
  dds <- DESeqDataSetFromHTSeqCount(
    sampleTable = sample_table,
    directory = count_dir,
    design = design
  )
  dds <- dds[rownames(dds) %in% gene_info$gene_id, ]
  dds <- dds[rowMeans(counts(dds)) >= min_count_mean, ]
  dds <- DESeq(dds)

  vsd <- vst(dds, blind = FALSE)
  vsd_tbl <- as.data.frame(assay(vsd)) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    right_join(gene_info, by = "gene_id") %>%
    relocate(gene_id, gene_name, gene_type)
  write_csv(vsd_tbl, file.path(outdir, "vsd.csv"))

  plot_deseq2_pca(
    vsd = vsd,
    output_pdf = file.path(outdir, "PCA_DESeq2.pdf"),
    output_csv = file.path(outdir, "PCA_DESeq2_data.csv"),
    intgroup = "group_name"
  )

  group_levels <- levels(colData(dds)$group_name)
  if (!(ref_group %in% group_levels)) {
    stop("ref_group is not in sample_table$group_name")
  }

  de_pair <- run_pairwise_de(dds, group_levels, ref_group, gene_info, outdir)
  de_custom <- run_custom_contrasts(dds, contrast_file, gene_info, outdir)
  de_result_list <- c(de_pair, de_custom)

  if (length(de_result_list) == 0) {
    stop("No DE comparison generated. Check group labels or contrast_file.")
  }

  saveRDS(de_result_list, file.path(outdir, "de_result_list.rds"))

  sig_de_genes <- lapply(
    de_result_list,
    function(x) {
      x %>%
        filter(!is.na(padj), abs(log2FoldChange) >= log2_cutoff, padj <= padj_cutoff) %>%
        pull(gene_id)
    }
  ) %>% unlist() %>% unique()

  write_lines(sig_de_genes, file.path(outdir, "significant_de_genes.txt"))

  if (run_volcano_plot) {
    run_volcano(de_result_list, padj_cutoff, log2_cutoff, outdir)
  }

  if (length(group_levels) > 2) {
    vsd_mat <- assay(vsd)
    run_timecourse_cluster(vsd_mat, sample_table, sig_de_genes, cluster_num, outdir)
  } else {
    message("Group number <= 2; skip time-course clustering.")
  }

  message("Analysis finished. Results in: ", outdir)
}

main()
