# RNAseqAnalysis

Question-driven RNA-seq differential expression analysis skill built around DESeq2, optional TCseq trend clustering, clusterProfiler enrichment, parameter-variant comparison, pre-run plan display, and markdown report generation.

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

- [references/skill_input.schema.json](/Users/xiaoclaw/DEanalysis/references/skill_input.schema.json)
- [references/skill_output.schema.json](/Users/xiaoclaw/DEanalysis/references/skill_output.schema.json)
- [references/intake_questionnaire.md](/Users/xiaoclaw/DEanalysis/references/intake_questionnaire.md)
- [references/report_template.md](/Users/xiaoclaw/DEanalysis/references/report_template.md)

## Output Layout

Each run writes:

- `resolved_config.json`
- `variant_summary.csv`
- `manifest.json`
- `RNAseqAnalysis_report.md`
- one subdirectory per parameter variant

Each variant directory contains:

- `vsd.csv`
- `PCA_DESeq2.pdf`
- `PCA_DESeq2_data.csv`
- `*_DE.csv`
- `*_summary.csv`
- optional `enrichment/`
- optional `tcseq/`

## Notes

- The current implementation prioritizes modular interfaces and extensibility.
- Enrichment uses GO BP and KEGG by default when the required annotation packages are available.
- TCseq is triggered automatically when group count exceeds the configured threshold or the study config explicitly requests trend analysis.
