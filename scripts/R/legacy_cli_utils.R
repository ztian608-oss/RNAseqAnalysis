print_help <- function() {
  cat(
    paste(
      "DESeq2 generic time-course pipeline",
      "",
      "Required arguments:",
      "  --count_dir=PATH          Directory containing HTSeq count files",
      "  --sample_table=FILE       CSV with sampleName, sampleFile, group_name",
      "  --gtf=FILE                GTF/GFF file with gene annotations",
      "  --ref_group=GROUP         Reference group for pairwise comparisons",
      "",
      "Optional arguments:",
      "  --outdir=DIR              Output directory (default: DESeq2_result)",
      "  --design='~ group_name'   Design formula (default: ~ group_name)",
      "  --contrast_file=FILE      CSV with columns treat,base for custom contrasts",
      "  --comparison_mode=MODE    pairwise | custom | both (default: both)",
      "  --padj_cutoff=NUM         Significant adjusted p-value (default: 0.01)",
      "  --log2_cutoff=NUM         Significant |log2FC| (default: 1)",
      "  --min_count_mean=NUM      Filter low-expression genes (default: 10)",
      "  --cluster_num=NUM         K for time-course kmeans (default: 6)",
      "  --seed=NUM                Random seed (default: 12345)",
      "  --run_volcano=TRUE/FALSE  Whether to draw volcano plots (default: TRUE)",
      "  --run_heatmap=TRUE/FALSE  Whether to draw DE heatmap (default: TRUE)",
      "  --help                    Show this help message",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    print_help()
    quit(save = "no", status = 0)
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

assert_file_exists <- function(path, desc) {
  if (!file.exists(path)) {
    stop(desc, " not found: ", path)
  }
}
