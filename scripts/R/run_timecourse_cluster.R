#!/usr/bin/env Rscript

aggregate_group_means <- function(vsd_mat, sample_table, ordered_groups) {
  group_means <- lapply(ordered_groups, function(group_label) {
    idx <- sample_table$sampleName[sample_table$group_name == group_label]
    rowMeans(vsd_mat[, idx, drop = FALSE], na.rm = TRUE)
  })
  pattern_mat <- do.call(cbind, group_means)
  colnames(pattern_mat) <- ordered_groups
  pattern_mat
}

safe_scale_rows <- function(mat) {
  if (nrow(mat) == 0) {
    return(mat)
  }
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  scaled <- sweep(mat, 1, mu, "-")
  scaled <- sweep(scaled, 1, sdv, "/")
  scaled
}

infer_cluster_trend <- function(centers) {
  first_point <- centers[, 1]
  last_point <- centers[, ncol(centers)]
  delta <- last_point - first_point
  dplyr::case_when(
    delta > 0.2 ~ "up",
    delta < -0.2 ~ "down",
    TRUE ~ "stable"
  )
}

reorder_tcseq_clusters <- function(tca, trend_table) {
  trend_order <- c("up", "stable", "down")
  ordered_clusters <- trend_table %>%
    dplyr::mutate(trend = factor(.data$trend, levels = trend_order)) %>%
    dplyr::arrange(.data$trend, dplyr::desc(.data$delta), .data$cluster) %>%
    dplyr::pull(.data$cluster) %>%
    as.character()

  mapping <- stats::setNames(seq_along(ordered_clusters), ordered_clusters)
  old_cluster <- as.character(tca@cluster)
  tca@cluster <- as.integer(mapping[old_cluster])
  tca@membership <- tca@membership[, ordered_clusters, drop = FALSE]
  colnames(tca@membership) <- as.character(seq_along(ordered_clusters))
  tca@centers <- tca@centers[ordered_clusters, , drop = FALSE]
  rownames(tca@centers) <- as.character(seq_along(ordered_clusters))
  tca
}

export_cluster_gene_lists <- function(cluster_assign, outdir) {
  dir.create(file.path(outdir, "timecourse_cluster_gene_lists"), showWarnings = FALSE, recursive = TRUE)
  split_df <- split(cluster_assign$gene_id, cluster_assign$cluster_label)
  for (nm in names(split_df)) {
    out_file <- file.path(outdir, "timecourse_cluster_gene_lists", paste0(nm, ".txt"))
    readr::write_lines(split_df[[nm]], out_file)
  }
}

run_timecourse_cluster <- function(vsd_mat, sample_table, sig_genes, ref_group, cluster_num, outdir, seed = 12345) {
  required_pkgs <- c("TCseq", "dplyr", "tidyr", "gridExtra", "readr", "ggplot2", "tibble")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages for time-course clustering: ", paste(missing_pkgs, collapse = ", "))
  }

  if (length(sig_genes) == 0) {
    message("No significant DE genes found; skip TCseq time-course clustering.")
    return(invisible(NULL))
  }

  if (!(ref_group %in% unique(as.character(sample_table$group_name)))) {
    stop("ref_group is not in sample_table$group_name for time-course clustering.")
  }

  sig_genes <- intersect(sig_genes, rownames(vsd_mat))
  if (length(sig_genes) < 2) {
    message("Fewer than 2 significant genes overlap with VST matrix; skip TCseq clustering.")
    return(invisible(NULL))
  }

  group_levels <- unique(as.character(sample_table$group_name))
  ordered_groups <- c(ref_group, setdiff(group_levels, ref_group))
  max_clusters <- max(2L, length(ordered_groups) * 2L)

  pattern_mat <- aggregate_group_means(vsd_mat[sig_genes, , drop = FALSE], sample_table, ordered_groups)
  pattern_scaled <- safe_scale_rows(pattern_mat)

  keep <- apply(pattern_scaled, 1, function(x) !any(is.na(x) | is.infinite(x)))
  pattern_scaled <- pattern_scaled[keep, , drop = FALSE]
  if (nrow(pattern_scaled) < 2) {
    message("Not enough genes after scaling/cleaning for TCseq clustering; skip.")
    return(invisible(NULL))
  }
  requested_cluster_num <- as.integer(cluster_num)
  cluster_num <- max(2L, min(requested_cluster_num, nrow(pattern_scaled), max_clusters))
  if (!is.na(requested_cluster_num) && requested_cluster_num > max_clusters) {
    message("cluster_num capped to ", cluster_num, " (<= 2x group count).")
  }

  set.seed(seed)
  tca <- TCseq::timeclust(pattern_scaled, algo = "cm", k = cluster_num)

  centers <- tca@centers
  trend <- infer_cluster_trend(centers)
  delta <- centers[, ncol(centers)] - centers[, 1]
  trend_table <- tibble::tibble(
    cluster = as.integer(seq_len(nrow(centers))),
    trend = trend,
    delta = as.numeric(delta)
  )
  tca <- reorder_tcseq_clusters(tca, trend_table)

  trend_table <- tibble::tibble(
    cluster = as.integer(seq_len(nrow(tca@centers))),
    trend = infer_cluster_trend(tca@centers),
    delta = as.numeric(tca@centers[, ncol(tca@centers)] - tca@centers[, 1])
  )

  grDevices::pdf(file.path(outdir, "timecourse_cluster_patterns.pdf"), width = 6, height = 6, onefile = TRUE)
  p <- TCseq::timeclustplot(
    tca,
    cols = 2,
    membership.color = c("royalblue", "moccasin", "red"),
    cl.color = "black",
    title.size = 8,
    axis.text.size = 8,
    legend.title.size = 8,
    legend.text.size = 8,
    axis.line.size = 0.4,
    axis.title.size = 8
  )
  gridExtra::grid.arrange(grobs = p, ncol = 2)
  grDevices::dev.off()

  cluster_assign <- tibble::tibble(
    gene_id = rownames(pattern_scaled),
    cluster = as.integer(tca@cluster)
  ) %>%
    dplyr::left_join(trend_table, by = c("cluster" = "cluster")) %>%
    dplyr::mutate(cluster_label = paste0("cluster_", .data$cluster, "_", .data$trend))

  readr::write_csv(cluster_assign, file.path(outdir, "timecourse_cluster_result.csv"))
  readr::write_csv(trend_table, file.path(outdir, "timecourse_cluster_trend_summary.csv"))
  export_cluster_gene_lists(cluster_assign, outdir)

  centers_long <- as.data.frame(tca@centers) %>%
    tibble::rownames_to_column("cluster") %>%
    dplyr::mutate(cluster = as.integer(.data$cluster)) %>%
    dplyr::left_join(trend_table, by = c("cluster" = "cluster")) %>%
    tidyr::pivot_longer(cols = all_of(ordered_groups), names_to = "group_name", values_to = "z_expr")

  p2 <- ggplot2::ggplot(
    centers_long,
    ggplot2::aes(x = .data$group_name, y = .data$z_expr, group = .data$cluster, color = .data$trend)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "TCseq cluster trend patterns", x = "Group", y = "Scaled expression") +
    ggplot2::facet_wrap(~cluster, scales = "free_y")

  ggplot2::ggsave(file.path(outdir, "timecourse_cluster_trend_lines.pdf"), plot = p2, width = 8, height = 5)

  invisible(list(
    tca = tca,
    trend_table = trend_table,
    cluster_assign = cluster_assign
  ))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  usage <- paste(
    "Usage:",
    "Rscript scripts/R/run_timecourse_cluster.R",
    "--vsd_csv <vsd.csv>",
    "--sample_table <sample_table.csv>",
    "--sig_genes <significant_de_genes.txt>",
    "--ref_group <control>",
    "--cluster_num <6>",
    "--outdir <output_dir>"
  )

  arg_value <- function(flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) == 0 || idx == length(args)) {
      return(default)
    }
    args[idx + 1]
  }

  required_pkgs <- c("readr", "dplyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages for CLI mode: ", paste(missing_pkgs, collapse = ", "))
  }

  vsd_csv <- arg_value("--vsd_csv")
  sample_table_file <- arg_value("--sample_table")
  sig_genes_file <- arg_value("--sig_genes")
  ref_group <- arg_value("--ref_group")
  cluster_num <- as.integer(arg_value("--cluster_num", "6"))
  outdir <- arg_value("--outdir", ".")
  seed <- as.integer(arg_value("--seed", "12345"))

  if (any(vapply(list(vsd_csv, sample_table_file, sig_genes_file, ref_group), is.null, logical(1)))) {
    stop(usage)
  }
  if (is.na(cluster_num) || cluster_num < 2) {
    stop("--cluster_num must be an integer >= 2")
  }

  vsd_tbl <- readr::read_csv(vsd_csv, show_col_types = FALSE)
  if (!("gene_id" %in% colnames(vsd_tbl))) {
    stop("vsd_csv must contain a 'gene_id' column.")
  }
  sample_cols <- setdiff(colnames(vsd_tbl), c("gene_id", "gene_name", "gene_type"))
  vsd_mat <- as.matrix(vsd_tbl[, sample_cols, drop = FALSE])
  rownames(vsd_mat) <- vsd_tbl$gene_id
  storage.mode(vsd_mat) <- "numeric"

  sample_table <- readr::read_csv(sample_table_file, show_col_types = FALSE) %>%
    dplyr::mutate(group_name = factor(.data$group_name))
  sig_genes <- readr::read_lines(sig_genes_file)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  run_timecourse_cluster(
    vsd_mat = vsd_mat,
    sample_table = sample_table,
    sig_genes = sig_genes,
    ref_group = ref_group,
    cluster_num = cluster_num,
    outdir = outdir,
    seed = seed
  )
}
