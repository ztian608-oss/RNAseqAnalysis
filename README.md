# DEanalysis

A modular DESeq2 pipeline for RNA-seq differential expression analysis with command-line execution, standardized outputs, and skill-ready templates.

---

## What this repository provides

- Differential expression analysis using **DESeq2** from HTSeq count files.
- Gene annotation mapping from **GTF/GFF** (`gene_id`, `gene_name`, optional `gene_type`).
- Flexible comparison modes:
  - pairwise contrasts versus a reference group,
  - custom contrasts from a `treat,base` CSV.
- Visualization and pattern analysis:
  - PCA,
  - volcano plots,
  - significant DE heatmap,
  - time-course/multi-group k-means clustering.
- Skill-oriented templates:
  - runnable command template,
  - JSON schemas for input/output contracts.

---

## Repository structure

- `scripts/deseq2_generic.R` — main entrypoint (orchestrates the full workflow).
- `scripts/R/utils_cli.R` — CLI parsing, argument helpers, and basic validations.
- `scripts/R/io.R` — sample table and GTF loading/validation.
- `scripts/R/analysis.R` — DESeq2 model construction, contrasts, result assembly, significant gene extraction.
- `scripts/R/plots.R` — PCA and volcano plotting.
- `scripts/R/clustering.R` — significant-gene heatmap and time-course clustering.
- `scripts/run_deseq2_generic.sh` — executable command template.
- `references/skill_input.schema.json` — skill input contract schema.
- `references/skill_output.schema.json` — skill output contract schema.
- `SKILL.md` + `agents/openai.yaml` — OpenClaw/Codex skill metadata.

---

## Input requirements

### Required CLI inputs

- `--count_dir`: directory containing HTSeq count files.
- `--sample_table`: CSV with required columns:
  - `sampleName`
  - `sampleFile`
  - `group_name`
- `--gtf`: annotation file (GTF/GFF).
- `--ref_group`: reference group name.

### Optional sample table column

- `group_order`: controls factor order for plotting/comparison order.

### Optional custom contrasts file

- `--contrast_file`: CSV with columns `treat,base`.
- Required when `--comparison_mode=custom` or `--comparison_mode=both` and custom contrasts are expected.

---

## Main CLI usage

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --outdir=DESeq2_result
```

Show help:

```bash
Rscript scripts/deseq2_generic.R --help
```

Use template (edit paths first):

```bash
bash scripts/run_deseq2_generic.sh
```

---

## Step-by-step runnable examples

### Step 1 — Base run (PCA + pairwise DE, no volcano/heatmap)

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --outdir=result_step1 \
  --comparison_mode=pairwise \
  --run_volcano=FALSE \
  --run_heatmap=FALSE
```

### Step 2 — Enable volcano and significant-gene heatmap

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --outdir=result_step2 \
  --comparison_mode=pairwise \
  --padj_cutoff=0.01 \
  --log2_cutoff=1 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE
```

### Step 3 — Add custom contrasts (optional)

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=control \
  --contrast_file=contrasts.csv \
  --comparison_mode=both \
  --outdir=result_step3 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE
```

### Step 4 — Multi-group / time-course clustering (group count > 2)

```bash
Rscript scripts/deseq2_generic.R \
  --count_dir=counts \
  --sample_table=sample_table.csv \
  --gtf=annotation.gtf \
  --ref_group=T0 \
  --outdir=result_step4 \
  --comparison_mode=both \
  --cluster_num=6 \
  --seed=12345 \
  --run_volcano=TRUE \
  --run_heatmap=TRUE
```

---

## Outputs

Core outputs:

- `vsd.csv`
- `PCA_DESeq2.pdf`
- `PCA_DESeq2_data.csv`
- `*_DE.csv`
- `*_summary.csv`
- `de_result_list.rds`
- `significant_de_genes.txt`

Optional/conditional outputs:

- `volcanoplot_*.pdf` (when `--run_volcano=TRUE`)
- `significant_de_heatmap.pdf` (when `--run_heatmap=TRUE`)
- `timecourse_cluster_heatmap.pdf`, `timecourse_cluster_result.csv`, `timecourse_cluster_patterns.pdf` (when group count > 2)

---

## PCA variance details

PCA output includes explicit variance metadata:

- Axis labels include explained variance for PC1 and PC2.
- Plot subtitle includes cumulative PC1+PC2 explained variance.
- `PCA_DESeq2_data.csv` includes:
  - `PC1_var_percent`
  - `PC2_var_percent`
  - `PC1_PC2_total_percent`

---

## Skill integration notes

For OpenClaw/agent workflow automation, use:

- `references/skill_input.schema.json` for input validation.
- `references/skill_output.schema.json` for output-port consistency.

This enables stable contracts between execution and downstream skill steps.

---

## Notes and guardrails

- The pipeline validates required files before analysis starts.
- If no genes remain after annotation/filtering, the run stops with an explicit error.
- If no significant genes are found, heatmap and clustering steps are skipped safely.
