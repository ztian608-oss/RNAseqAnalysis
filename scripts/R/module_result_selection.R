select_final_variant <- function(stability_tbl) {
  stability_tbl |>
    dplyr::arrange(
      dplyr::desc(score),
      dplyr::desc(comparison_overlap_mean),
      dplyr::desc(marker_direction_consistency),
      dplyr::desc(top_gene_jaccard),
      total_sig
    ) |>
    dplyr::slice(1)
}
