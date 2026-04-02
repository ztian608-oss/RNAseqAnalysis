---
name: rnaseq-analysis
description: Plan and run modular RNA-seq differential expression workflows for bulk count data. Use when users need DESeq2-based analysis driven by research questions, including pre-run analysis planning, optional background-question prompts, sample grouping priorities, threshold preferences, optional TCseq trend clustering for more than 3 groups, clusterProfiler enrichment, multi-parameter stability comparison, and markdown report generation with standardized outputs.
---

# RNAseqAnalysis

Use this skill when the user wants an end-to-end differential expression workflow that adapts to the biological question instead of only running a fixed CLI.

This skill is intended to behave like a professional bioinformatics analysis service, not a one-off script. The priorities are:

1. interpretation quality,
2. parameter discipline,
3. reproducibility,
4. extensibility for planner or agent orchestration.

## Workflow

1. Before execution, show the planned analysis steps and expected outputs.
2. Collect or infer the user context before analysis. Ask for:
   - research objective,
   - disease, pathway, or marker focus,
   - target sample groups and preferred contrasts,
   - threshold preferences,
   - whether trend analysis is needed,
   - whether to favor statistical strictness or biological interpretability.
3. If the prompt contains no usable biological background, ask the user for study background before selecting the final result variant.
4. Convert the answers into a normalized config. Use the schema in `references/skill_input.schema.json`.
5. Validate user inputs:
   - count matrix or HTSeq count directory,
   - metadata table,
   - optional contrast file,
   - optional annotation and background gene sets.
   - parameter policy in `references/parameter_policy.md`
6. Run DESeq2 analysis through `scripts/run_rnaseq_analysis.R`.
7. If the effective group count is greater than 3, or the config explicitly enables trend analysis, run the TCseq module.
8. Run clusterProfiler enrichment on selected DEG sets and on TCseq clusters when available.
9. If multiple threshold or clustering parameter sets are configured, compare stability and select the final result bundle.
10. Generate a markdown report and a machine-readable output manifest.

## Key Files

- `scripts/run_rnaseq_analysis.R`: orchestrates the modular workflow.
- `scripts/R/module_preflight.R`: pre-run plan display and missing-background prompts.
- `scripts/R/module_config.R`: config loading and variant expansion.
- `scripts/R/module_intake.R`: intake questions and user-context normalization.
- `scripts/R/module_io.R`: input loading.
- `scripts/R/module_validation.R`: schema-like validation and initialization.
- `scripts/R/module_design.R`: design and contrast planning.
- `scripts/R/module_deseq2_engine.R`: DESeq2 object construction and fitting.
- `scripts/R/module_deg_analysis.R`: results extraction, shrinkage, DEG filtering.
- `scripts/R/module_biological_interpretation.R`: rule-based biological explanation layer.
- `scripts/R/module_tcseq_analysis.R`: trend clustering.
- `scripts/R/module_enrichment.R`: clusterProfiler-based ORA/GSEA wrappers.
- `scripts/R/module_stability.R`: multi-parameter comparison and scoring.
- `scripts/R/module_report_generation.R`: markdown report rendering.
- `scripts/R/module_workflow_orchestrator.R`: workflow execution and manifest export.

## References

Before any real run, read `references/intake_questionnaire.md` and check whether the user already provided:
- disease or study context,
- pathway or marker priorities,
- why a particular variant should be favored.

- Read `references/intake_questionnaire.md` before collecting user requirements.
- Read `references/parameter_policy.md` before accepting any threshold or clustering parameters.
- Read `references/skill_architecture.md` when modifying workflow boundaries or planner integration.
- Read `references/interpretation_rules.json` when extending biological explanation panels or context tags.
- Read `references/report_template.md` when modifying report structure.
- Use `references/skill_input.schema.json` and `references/skill_output.schema.json` for interface contracts.

## Execution Pattern

Primary run:

```bash
Rscript scripts/run_rnaseq_analysis.R --config=templates/rnaseq_analysis_full_config.json
```

Compatibility run for the legacy DESeq2-only CLI:

```bash
Rscript scripts/legacy_deseq2_cli.R --count_dir=counts --sample_table=samples.csv --gtf=genes.gtf --ref_group=control --outdir=result
```

## Guardrails

- Prefer `DESeqDataSetFromMatrix` when the user already has a count matrix.
- Prefer `DESeqDataSetFromHTSeqCount` only for HTSeq-style per-sample files.
- Do not run time-course clustering for 2-group comparisons unless the user explicitly requests trend analysis.
- Treat enrichment as question-dependent interpretation, not an automatic substitute for effect-size review.
- When multiple parameter variants are run, report both the chosen final result and the competing variants.
- Only accept `log2FC` thresholds from `{1, 1.5, 2}`.
- Only accept `padj` thresholds from `{0.05, 0.01}`.
- Reject ad hoc “soft” thresholds even if the user asks for them, unless the user explicitly wants a non-standard exploratory run and the report labels it as such.
