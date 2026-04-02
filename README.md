# RNAseqAnalysis

Question-driven RNA-seq differential expression analysis skill built around DESeq2, optional TCseq trend clustering, clusterProfiler enrichment, parameter-variant comparison, pre-run plan display, and markdown report generation.

The current workflow applies a low-expression prefilter before DESeq2 fitting, interprets differential expression relative to the configured reference group, and reports direction-specific GO/KEGG enrichment for up-regulated and down-regulated genes separately.

## Current Design

- `scripts/run_rnaseq_analysis.R` is the main config-driven entrypoint.
- `scripts/legacy_deseq2_cli.R` is a compatibility wrapper for the previous CLI.
- `scripts/R/` contains modular workflow components:
  - preflight planning and background prompts,
  - config and intake normalization,
  - input loading and validation,
  - design and contrast planning,
  - DESeq2 fitting and DEG extraction,
  - TCseq trend clustering,
  - enrichment,
  - stability scoring,
  - report and manifest export.

## Main Run

```bash
Rscript scripts/run_rnaseq_analysis.R --config=templates/rnaseq_analysis_full_config.json
```

## Legacy-Compatible Run

```bash
Rscript scripts/legacy_deseq2_cli.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --outdir=RNAseqAnalysis_output
```

## Config Model

The normalized config contains these sections:

- `study`
- `input`
- `design`
- `thresholds`
- `shrinkage`
- `trend`
- `enrichment`
- `output`

See:

- [references/skill_input.schema.json](/Users/xiaoclaw/RNAseqAnalysis/references/skill_input.schema.json)
- [references/skill_output.schema.json](/Users/xiaoclaw/RNAseqAnalysis/references/skill_output.schema.json)
- [references/intake_questionnaire.md](/Users/xiaoclaw/RNAseqAnalysis/references/intake_questionnaire.md)
- [references/report_template.md](/Users/xiaoclaw/RNAseqAnalysis/references/report_template.md)

## Output Layout

Each run writes:

- `config/resolved_config.json`
- `comparison/variant_summary.csv`
- `manifest.json`
- `report/RNAseqAnalysis_report.md`
- variant-scoped outputs under `deseq2/`, `figures/`, `tables/`, `comparison/`, `enrichment/`, and `tcseq/`

For the selected variant, the report links to files such as:

- `deseq2/<variant>/vsd.csv`
- `figures/<variant>/PCA_DESeq2.pdf`
- `tables/<variant>/PCA_DESeq2_data.csv`
- `tables/<variant>/*_DE.csv`
- `comparison/<variant>/*_summary.csv`
- `enrichment/<variant>/<comparison>_up_go_bp.csv`
- `enrichment/<variant>/<comparison>_up_kegg.csv`
- `enrichment/<variant>/<comparison>_down_go_bp.csv`
- `enrichment/<variant>/<comparison>_down_kegg.csv`
- `*_DE.csv`
- `*_summary.csv`
- optional `tcseq/<variant>/`

## Notes

- The current implementation prioritizes modular interfaces and extensibility.
- Enrichment uses GO BP and KEGG by default when the required annotation packages are available.
- `min_count_mean` defaults to `10` and must be greater than `0`.
- TCseq is triggered automatically when group count exceeds the configured threshold or the study config explicitly requests trend analysis.
