scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  s[s == 0] <- 1
  (x - m) / s
}

build_dds <- function(sample_table, count_dir, design_formula, gene_info, min_count_mean) {
  design <- as.formula(design_formula)

  dds <- DESeq2::DESeqDataSetFromHTSeqCount(
    sampleTable = sample_table,
    directory = count_dir,
    design = design
  )

  dds <- dds[rownames(dds) %in% gene_info$gene_id, ]
  dds <- dds[rowMeans(DESeq2::counts(dds)) >= min_count_mean, ]

  if (nrow(dds) == 0) {
    stop("No genes remain after gene_id mapping + min_count_mean filter.")
  }

  DESeq2::DESeq(dds)
}

assemble_result_table <- function(vsd_sub, res, gene_info) {
  cbind(as.data.frame(SummarizedExperiment::assay(vsd_sub)), as.data.frame(res)) %>%
    tibble::rownames_to_column("gene_id") %>%
    tibble::as_tibble() %>%
    dplyr::right_join(gene_info, by = "gene_id") %>%
    dplyr::relocate(gene_id, gene_name, gene_type)
}

run_single_contrast <- function(dds, treat, base, gene_info, outdir) {
  dds_sub <- dds[, dds$group_name %in% c(base, treat)]
  dds_sub$group_name <- droplevels(dds_sub$group_name)

  if (nlevels(dds_sub$group_name) < 2) {
    warning("Skip contrast ", treat, " vs ", base, ": insufficient group levels")
    return(NULL)
  }

  dds_sub <- DESeq2::DESeq(dds_sub)
  res <- DESeq2::results(dds_sub, contrast = c("group_name", treat, base))
  vsd_sub <- DESeq2::vst(dds_sub, blind = FALSE)

  res_tbl <- assemble_result_table(vsd_sub, res, gene_info)
  nm <- paste0(treat, "_vs_", base)
  readr::write_csv(res_tbl, file.path(outdir, paste0(nm, "_DE.csv")))

  sig_tbl <- res_tbl %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::summarise(
      total = dplyr::n(),
      up = sum(log2FoldChange > 0 & padj <= 0.05),
      down = sum(log2FoldChange < 0 & padj <= 0.05)
    ) %>%
    dplyr::mutate(comparison = nm)

  readr::write_csv(sig_tbl, file.path(outdir, paste0(nm, "_summary.csv")))
  list(name = nm, table = res_tbl)
}

run_pairwise_de <- function(dds, group_levels, ref_group, gene_info, outdir) {
  out <- list()
  for (cmp in group_levels[group_levels != ref_group]) {
    one <- run_single_contrast(dds, cmp, ref_group, gene_info, outdir)
    if (!is.null(one)) out[[one$name]] <- one$table
  }
  out
}

run_custom_contrasts <- function(dds, contrast_file, gene_info, outdir) {
  if (is.null(contrast_file) || contrast_file == "") return(list())

  assert_file_exists(contrast_file, "contrast_file")
  contrast_tbl <- readr::read_csv(contrast_file, show_col_types = FALSE)
  required_cols <- c("treat", "base")
  if (!all(required_cols %in% colnames(contrast_tbl))) {
    stop("contrast_file must contain columns: treat, base")
  }

  out <- list()
  for (i in seq_len(nrow(contrast_tbl))) {
    one <- run_single_contrast(
      dds = dds,
      treat = contrast_tbl$treat[[i]],
      base = contrast_tbl$base[[i]],
      gene_info = gene_info,
      outdir = outdir
    )
    if (!is.null(one)) out[[one$name]] <- one$table
  }

  out
}

extract_significant_genes <- function(de_result_list, log2_cutoff, padj_cutoff) {
  lapply(
    de_result_list,
    function(x) {
      x %>%
        dplyr::filter(!is.na(padj), abs(log2FoldChange) >= log2_cutoff, padj <= padj_cutoff) %>%
        dplyr::pull(gene_id)
    }
  ) %>% unlist() %>% unique()
}
