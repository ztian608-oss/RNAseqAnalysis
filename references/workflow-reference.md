# Workflow Reference

## Required sample table schema

Required columns:
- `sampleName`
- `sampleFile`
- `group_name`

Optional:
- `group_order`

## Contrast file schema

CSV columns:
- `treat`
- `base`

## Important runtime rules

- Use `comparison_mode=pairwise` when no custom contrast file is provided.
- Use `comparison_mode=custom` or `both` only with a valid `contrast_file`.
- Time-course clustering is only meaningful when group count is greater than 2.

## Key output artifacts

Core:
- `vsd.csv`
- `PCA_DESeq2.pdf`
- `PCA_DESeq2_data.csv`
- `manifest.json`
- `RNAseqAnalysis_report.md`
- `variant_summary.csv`

Per-comparison:
- `*_DE.csv`
- `*_summary.csv`
- `*_up_go_bp.csv`
- `*_up_kegg.csv`
- `*_down_go_bp.csv`
- `*_down_kegg.csv`

Optional:
- `volcano_*.pdf`
- `significant_de_heatmap.pdf`
- `timecourse_cluster_heatmap.pdf`
- `timecourse_cluster_result.csv`
- `timecourse_cluster_patterns.pdf`
