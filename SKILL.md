---
name: deseq2-generic-analysis
description: Execute and adapt the modular DESeq2 RNA-seq differential expression workflow in this repository. Use when users need to run DESeq2 from HTSeq counts, configure pairwise/custom contrasts, generate PCA/volcano/heatmap/time-course clustering outputs, or package the workflow into reproducible skill/API-style runs with stable JSON input/output contracts.
---

# DESeq2 Generic Analysis Skill

Run the repository pipeline with correct inputs, robust parameter choices, and standardized output reporting.

## Workflow

1. Validate required inputs:
   - count directory,
   - sample table (`sampleName`, `sampleFile`, `group_name`),
   - GTF/GFF,
   - reference group.
2. Choose comparison mode:
   - `pairwise`: compare each non-reference group vs reference,
   - `custom`: only contrast file comparisons,
   - `both`: pairwise + custom.
3. Run `scripts/deseq2_generic.R` with selected options.
4. Summarize key outputs for downstream steps.
5. When requested, emit interface-aligned JSON using schema templates in `references/`.

## Command patterns

### Base run

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=<counts_dir> \
  --sample_table=<sample_table.csv> \
  --gtf=<annotation.gtf> \
  --ref_group=<reference_group> \
  --outdir=<output_dir> \
  --comparison_mode=pairwise
```

### Run with custom contrasts

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=<counts_dir> \
  --sample_table=<sample_table.csv> \
  --gtf=<annotation.gtf> \
  --ref_group=<reference_group> \
  --contrast_file=<contrasts.csv> \
  --comparison_mode=both \
  --outdir=<output_dir>
```

## Output checklist

Always verify and report:

- `vsd.csv`
- `PCA_DESeq2.pdf`
- `PCA_DESeq2_data.csv`
- `de_result_list.rds`
- `significant_de_genes.txt`
- any `*_DE.csv` and `*_summary.csv`
- optional artifacts (`volcanoplot_*.pdf`, `significant_de_heatmap.pdf`, time-course cluster files)

## References

- For CLI and artifact details, read `references/workflow-reference.md`.
- For strict input/output contracts, read:
  - `references/skill_input.schema.json`
  - `references/skill_output.schema.json`
