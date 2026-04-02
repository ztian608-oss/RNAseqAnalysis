#!/usr/bin/env bash
set -euo pipefail

# ===== Example command template for DESeq2 generic pipeline =====
# Copy this file and replace paths/groups for your project.

Rscript legacy_deseq2_cli.R \
  --count_dir=/path/to/htseq_counts \
  --sample_table=/path/to/sample_table.csv \
  --gtf=/path/to/annotation.gtf \
  --ref_group=control \
  --outdir=/path/to/output/RNAseqAnalysis_output \
  --design='~ group_name' \
  --comparison_mode=pairwise \
  --padj_cutoff=0.01 \
  --log2_cutoff=1 \
  --min_count_mean=10 \
  --cluster_num=6 \
  --seed=12345 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE

# If you need custom contrasts:
# 1) prepare contrasts.csv with columns: treat,base
# 2) set --comparison_mode=both (or custom)
# 3) add: --contrast_file=/path/to/contrasts.csv
