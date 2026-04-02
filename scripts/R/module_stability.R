jaccard_index <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

mean_or_default <- function(x, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  mean(x)
}

marker_direction_profile <- function(focus_marker_table) {
  if (is.null(focus_marker_table) || nrow(focus_marker_table) == 0) {
    return(character())
  }
  out <- apply(focus_marker_table, 1, function(x) {
    marker <- ifelse(is.na(x[["gene_name"]]) || x[["gene_name"]] == "", x[["gene_id"]], x[["gene_name"]])
    direction <- ifelse(as.numeric(x[["log2FoldChange"]]) >= 0, "up", "down")
    paste(x[["comparison"]], marker, direction, sep = "::")
  })
  unique(out)
}

comparison_overlap_score <- function(result_i, result_j) {
  comps <- intersect(names(result_i$significant_gene_sets), names(result_j$significant_gene_sets))
  if (length(comps) == 0) return(0)
  mean(vapply(comps, function(nm) {
    jaccard_index(result_i$significant_gene_sets[[nm]], result_j$significant_gene_sets[[nm]])
  }, numeric(1)))
}

clamp01 <- function(x) {
  pmax(0, pmin(1, x))
}

score_variant <- function(result, study_cfg, median_sig_count) {
  denom <- max(median_sig_count, 1)
  sig_count_stability <- clamp01(1 - abs(result$total_sig - median_sig_count) / denom)

  expected_marker_hits <- length(study_cfg$focus_markers %||% character()) * length(result$comparisons %||% character())
  focus_marker_coverage <- if (expected_marker_hits <= 0) {
    0
  } else {
    clamp01(result$focus_marker_hits / expected_marker_hits)
  }

  weights <- switch(
    study_cfg$emphasis,
    strict_stats = c(
      comparison_overlap_mean = 0.35,
      top_gene_jaccard = 0.20,
      marker_direction_consistency = 0.10,
      focus_marker_coverage = 0.10,
      sig_count_stability = 0.25
    ),
    biological_interpretation = c(
      comparison_overlap_mean = 0.20,
      top_gene_jaccard = 0.15,
      marker_direction_consistency = 0.25,
      focus_marker_coverage = 0.25,
      sig_count_stability = 0.15
    ),
    c(
      comparison_overlap_mean = 0.30,
      top_gene_jaccard = 0.20,
      marker_direction_consistency = 0.20,
      focus_marker_coverage = 0.15,
      sig_count_stability = 0.15
    )
  )

  list(
    score =
      weights[["comparison_overlap_mean"]] * result$comparison_overlap_mean +
      weights[["top_gene_jaccard"]] * result$top_gene_jaccard +
      weights[["marker_direction_consistency"]] * result$marker_direction_consistency +
      weights[["focus_marker_coverage"]] * focus_marker_coverage +
      weights[["sig_count_stability"]] * sig_count_stability,
    sig_count_stability = sig_count_stability,
    focus_marker_coverage = focus_marker_coverage
  )
}

compare_variants <- function(variant_results, study_cfg) {
  all_counts <- vapply(variant_results, function(x) x$total_sig, numeric(1))
  top_lists <- lapply(variant_results, function(x) x$top_genes)
  marker_profiles <- lapply(variant_results, function(x) marker_direction_profile(x$focus_marker_table))
  median_sig_count <- stats::median(all_counts)

  stability <- lapply(seq_along(variant_results), function(i) {
    others <- setdiff(seq_along(variant_results), i)
    top_overlap <- if (length(others) == 0) 1 else mean(vapply(others, function(j) {
      jaccard_index(top_lists[[i]], top_lists[[j]])
    }, numeric(1)))
    comp_overlap <- if (length(others) == 0) 1 else mean(vapply(others, function(j) {
      comparison_overlap_score(variant_results[[i]], variant_results[[j]])
    }, numeric(1)))
    marker_consistency <- if (length(others) == 0) 1 else mean(vapply(others, function(j) {
      jaccard_index(marker_profiles[[i]], marker_profiles[[j]])
    }, numeric(1)))
    focus_hits <- if (nrow(variant_results[[i]]$focus_marker_table) == 0) 0 else nrow(variant_results[[i]]$focus_marker_table)

    enriched <- within(variant_results[[i]], {
      top_gene_jaccard <- top_overlap
      comparison_overlap_mean <- comp_overlap
      marker_direction_consistency <- marker_consistency
      focus_marker_hits <- focus_hits
    })

    score_components <- score_variant(enriched, study_cfg, median_sig_count)

    list(
      variant_id = enriched$variant_id,
      total_sig = enriched$total_sig,
      top_gene_jaccard = enriched$top_gene_jaccard,
      comparison_overlap_mean = enriched$comparison_overlap_mean,
      marker_direction_consistency = enriched$marker_direction_consistency,
      focus_marker_hits = enriched$focus_marker_hits,
      focus_marker_coverage = score_components$focus_marker_coverage,
      sig_count_stability = score_components$sig_count_stability,
      score = score_components$score
    )
  })

  dplyr::bind_rows(stability)
}
