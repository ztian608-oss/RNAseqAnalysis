scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  s[s == 0 | is.na(s)] <- 1
  (x - m) / s
}

plot_deseq2_pca <- function(vsd, output_pdf, output_csv, intgroup = "group_name") {
  pca_data <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"), 2)
  pc12_total <- round(sum(percent_var[1:2]), 2)

  pca_data$PC1_var_percent <- percent_var[1]
  pca_data$PC2_var_percent <- percent_var[2]
  pca_data$PC1_PC2_total_percent <- pc12_total

  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(PC1, PC2, color = .data[[intgroup]], label = name)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::labs(
      title = "PCA (DESeq2 VST)",
      subtitle = paste0("PC1+PC2 explained variance: ", pc12_total, "%"),
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = intgroup,
      caption = "Variance values are derived from DESeq2::plotPCA percentVar"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(output_pdf, plot = p, width = 5, height = 4)
  readr::write_csv(pca_data, output_csv)
  pca_data
}

comparison_label_from_name <- function(nm, contrast = NULL) {
  # Keep DESeq2 direction consistency: "treatment vs reference"
  # For standard names (A_vs_B), A is the treatment group.
  if (!is.null(contrast) && length(contrast) >= 3) {
    return(paste0(as.character(contrast[[2]]), " vs ", as.character(contrast[[3]])))
  }

  if (grepl("_vs_", nm, fixed = TRUE)) {
    parts <- strsplit(nm, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      return(paste0(parts[[1]], " vs ", parts[[2]]))
    }
  }
  # Fallback for non-standard names, preserving original direction text where possible.
  if (grepl(" vs ", nm, fixed = TRUE)) {
    parts <- strsplit(nm, " vs ", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      return(paste0(parts[[1]], " vs ", parts[[2]]))
    }
  }
  nm
}

plot_volcano <- function(file_path, comparison_name, p_cutoff = 0.05, log2fc_cutoff = 0.1) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)

  needed_cols <- c("padj", "log2FoldChange")
  if (!all(needed_cols %in% colnames(df))) {
    stop("Volcano input is missing required columns: padj, log2FoldChange")
  }

  if (!"Gene" %in% colnames(df)) {
    if ("gene_name" %in% colnames(df)) {
      df$Gene <- df$gene_name
    } else if ("gene_id" %in% colnames(df)) {
      df$Gene <- df$gene_id
    } else {
      df$Gene <- rownames(df)
    }
  }

  df_plot <- df |>
    dplyr::mutate(
      min_nonzero_padj = suppressWarnings(min(padj[padj > 0], na.rm = TRUE)),
      padj_adj = dplyr::if_else(
        is.na(padj),
        NA_real_,
        dplyr::if_else(
          padj == 0,
          dplyr::if_else(is.finite(min_nonzero_padj), min_nonzero_padj * 0.1, 1e-300),
          padj
        )
      ),
      log_p = -log10(padj_adj),
      sig = dplyr::case_when(
        !is.na(padj) & padj < p_cutoff & log2FoldChange >= log2fc_cutoff ~ "up",
        !is.na(padj) & padj < p_cutoff & log2FoldChange <= -log2fc_cutoff ~ "down",
        TRUE ~ "non-sig"
      )
    ) |>
    dplyr::select(-min_nonzero_padj) |>
    dplyr::filter(!is.na(log2FoldChange), !is.na(log_p), is.finite(log_p))

  if (nrow(df_plot) == 0) {
    return(NULL)
  }

  max_fc <- max(abs(df_plot$log2FoldChange), na.rm = TRUE)

  up_top5 <- df_plot |>
    dplyr::filter(sig == "up") |>
    dplyr::slice_max(log2FoldChange, n = 5, with_ties = FALSE) |>
    dplyr::pull(Gene)
  down_top5 <- df_plot |>
    dplyr::filter(sig == "down") |>
    dplyr::slice_min(log2FoldChange, n = 5, with_ties = FALSE) |>
    dplyr::pull(Gene)
  top_genes <- unique(c(up_top5, down_top5))

  point_aes <- ggplot2::aes(size = baseMean)
  if (!"baseMean" %in% colnames(df_plot) || all(is.na(df_plot$baseMean))) {
    df_plot$baseMean <- 1
    point_aes <- ggplot2::aes()
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = log2FoldChange, y = log_p, color = sig)) +
    ggplot2::geom_point(point_aes, alpha = 0.65) +
    ggplot2::scale_color_manual(values = c("up" = "#C6295C", "down" = "#2C6DB2", "non-sig" = "grey80")) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::coord_cartesian(xlim = c(-max_fc, max_fc)) +
    ggplot2::labs(
      x = "log2(Fold Change)",
      y = "-log10(Adjusted P.Value)",
      title = paste0("Volcano Plot: ", comparison_name),
      color = "Significance"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if ("baseMean" %in% colnames(df_plot) && any(df_plot$baseMean > 0, na.rm = TRUE)) {
    p <- p + ggplot2::scale_size_continuous(range = c(0.5, 3), trans = "log10", guide = "none")
  }

  label_data <- df_plot |> dplyr::filter(Gene %in% top_genes)
  if (nrow(label_data) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = Gene),
        size = 3,
        fontface = "italic",
        segment.color = "black",
        max.overlaps = 20,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(label = Gene),
        size = 2.8,
        fontface = "italic",
        check_overlap = TRUE,
        show.legend = FALSE
      )
    }
  }

  p
}

run_volcano <- function(tbl, nm, padj_cutoff, log2_cutoff, outdir) {
  tbl <- tbl[!is.na(tbl$log2FoldChange), , drop = FALSE]
  if (nrow(tbl) == 0) return(NULL)

  out_file <- file.path(outdir, paste0("volcano_", nm, ".pdf"))
  tmp_csv <- tempfile(pattern = paste0("volcano_", nm, "_"), fileext = ".csv")

  volcano_tbl <- as.data.frame(tbl)
  if (!"Gene" %in% colnames(volcano_tbl)) {
    gene_name <- if ("gene_name" %in% colnames(volcano_tbl)) volcano_tbl$gene_name else NA_character_
    gene_id <- if ("gene_id" %in% colnames(volcano_tbl)) volcano_tbl$gene_id else rownames(volcano_tbl)
    volcano_tbl$Gene <- ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  }

  utils::write.csv(volcano_tbl, tmp_csv, row.names = FALSE)
  comparison_name <- comparison_label_from_name(nm, contrast = attr(tbl, "contrast"))
  p <- plot_volcano(tmp_csv, comparison_name, p_cutoff = padj_cutoff, log2fc_cutoff = log2_cutoff)
  if (is.null(p)) return(NULL)

  ggplot2::ggsave(out_file, plot = p, width = 6, height = 5)
  out_file
}

run_ma_plot <- function(tbl, nm, padj_cutoff, log2_cutoff, outdir) {
  if (!all(c("baseMean", "log2FoldChange", "padj") %in% colnames(tbl))) {
    return(NULL)
  }
  x <- tbl[!is.na(tbl$baseMean) & !is.na(tbl$log2FoldChange), , drop = FALSE]
  if (nrow(x) == 0) return(NULL)

  x$significance <- dplyr::case_when(
    !is.na(x$padj) & x$padj <= padj_cutoff & x$log2FoldChange >= log2_cutoff ~ "Up",
    !is.na(x$padj) & x$padj <= padj_cutoff & x$log2FoldChange <= -log2_cutoff ~ "Down",
    TRUE ~ "NS"
  )

  out_file <- file.path(outdir, paste0("MA_", nm, ".pdf"))
  p <- ggplot2::ggplot(x, ggplot2::aes(log10(baseMean + 1), log2FoldChange, color = significance)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.65) +
    ggplot2::geom_hline(yintercept = c(-log2_cutoff, log2_cutoff), linetype = 2, color = "grey50") +
    ggplot2::scale_color_manual(values = c(Up = "#b2182b", Down = "#2166ac", NS = "grey75")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("MA plot:", nm),
      x = "log10(baseMean + 1)",
      y = "log2 fold change"
    )
  ggplot2::ggsave(out_file, plot = p, width = 5.5, height = 5)
  out_file
}

run_sig_heatmap <- function(vsd_mat, sig_genes, metadata, sample_column, group_column, outdir) {
  sig_mat <- vsd_mat[intersect(sig_genes, rownames(vsd_mat)), , drop = FALSE]
  if (nrow(sig_mat) < 2) return(NULL)
  if (!requireNamespace("pheatmap", quietly = TRUE)) return(NULL)
  ann_col <- metadata[, c(sample_column, group_column), drop = FALSE]
  rownames(ann_col) <- ann_col[[sample_column]]
  ann_col[[sample_column]] <- NULL
  out_file <- file.path(outdir, "significant_de_heatmap.pdf")
  pheatmap::pheatmap(
    scale_rows(sig_mat),
    show_rownames = FALSE,
    annotation_col = ann_col,
    filename = out_file,
    width = 7,
    height = 8
  )
  out_file
}

plot_enrichment_bubble <- function(file, output_pdf, title = "GO Enrichment", top_n = 30, wrap_width = 45) {
  if (!file.exists(file)) return(NULL)
  df <- read.csv(file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)

  if (!"Description" %in% colnames(df)) return(NULL)

  p_col <- dplyr::case_when(
    "p.adjust" %in% colnames(df) ~ "p.adjust",
    "qvalue" %in% colnames(df) ~ "qvalue",
    "pvalue" %in% colnames(df) ~ "pvalue",
    TRUE ~ ""
  )
  if (!nzchar(p_col)) return(NULL)

  if (!"FoldEnrichment" %in% colnames(df)) {
    if (all(c("GeneRatio", "BgRatio") %in% colnames(df))) {
      parse_ratio <- function(x) {
        sapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
          if (length(z) != 2) return(NA_real_)
          as.numeric(z[[1]]) / as.numeric(z[[2]])
        })
      }
      df$FoldEnrichment <- parse_ratio(df$GeneRatio) / parse_ratio(df$BgRatio)
    } else {
      return(NULL)
    }
  }

  if (!"Count" %in% colnames(df)) {
    if ("geneID" %in% colnames(df)) {
      df$Count <- lengths(strsplit(as.character(df$geneID), "/", fixed = TRUE))
    } else {
      df$Count <- 1
    }
  }

  top_terms <- df |>
    dplyr::mutate(
      p_numeric = as.numeric(.data[[p_col]]),
      log_p = -log10(dplyr::if_else(is.na(p_numeric) | p_numeric <= 0, 1e-300, p_numeric)),
      FoldEnrichment = as.numeric(FoldEnrichment),
      Count = as.numeric(Count)
    ) |>
    dplyr::filter(is.finite(log_p), is.finite(FoldEnrichment), is.finite(Count)) |>
    dplyr::arrange(dplyr::desc(log_p)) |>
    dplyr::slice_head(n = top_n)

  if (nrow(top_terms) == 0) return(NULL)

  p <- ggplot2::ggplot(
    top_terms,
    ggplot2::aes(
      x = FoldEnrichment,
      y = stats::reorder(Description, FoldEnrichment)
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = log_p), alpha = 0.8) +
    ggplot2::scale_color_gradient(low = "#377EB8", high = "#E41A1C") +
    ggplot2::scale_size_continuous(
      range = c(3, 10),
      name = "Gene Count",
      breaks = scales::breaks_pretty(n = 4)
    ) +
    ggplot2::labs(
      title = title,
      x = "Fold Enrichment",
      y = NULL,
      color = expression(-log[10](pvalue))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(color = "black", size = 12),
      axis.text.y = ggplot2::element_text(color = "black", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30", hjust = 0.5),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(20, 20, 20, 20),
      legend.key.size = grid::unit(0.8, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(
        barwidth = grid::unit(0.5, "cm"),
        barheight = grid::unit(2, "cm")
      ),
      size = ggplot2::guide_legend(
        keywidth = grid::unit(0.5, "cm"),
        keyheight = grid::unit(0.2, "cm")
      )
    )

  ggplot2::ggsave(output_pdf, plot = p, width = 8, height = 6)
  output_pdf
}

run_enrichment_bubble_plots <- function(enrich_files, comparison_name, outdir) {
  if (is.null(enrich_files) || length(enrich_files) == 0) return(list())
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out <- list()
  label <- gsub("[^A-Za-z0-9_]+", "_", comparison_name)

  for (direction in c("up", "down")) {
    entry <- enrich_files[[direction]]
    if (is.null(entry)) next
    if (!is.null(entry$go_bp)) {
      out[[paste0(direction, "_go_bp")]] <- plot_enrichment_bubble(
        file = entry$go_bp,
        output_pdf = file.path(outdir, paste0("enrichment_bubble_", label, "_", direction, "_go_bp.pdf")),
        title = paste0(comparison_name, ": GO ", ifelse(direction == "up", "Up-regulated", "Down-regulated"))
      )
    }
    if (!is.null(entry$kegg)) {
      out[[paste0(direction, "_kegg")]] <- plot_enrichment_bubble(
        file = entry$kegg,
        output_pdf = file.path(outdir, paste0("enrichment_bubble_", label, "_", direction, "_kegg.pdf")),
        title = paste0(comparison_name, ": KEGG ", ifelse(direction == "up", "Up-regulated", "Down-regulated"))
      )
    }
  }
  out[!vapply(out, is.null, logical(1))]
}
