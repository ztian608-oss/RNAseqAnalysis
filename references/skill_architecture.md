# Skill Architecture

## Design Goal

RNAseqAnalysis is a reusable analysis skill, not a monolithic script. It should be easy to run from a human-authored config, but also easy for a planner or agent system to call in stages.

## Layered Design

### 1. Intake Layer

Responsible for:

- collecting study background
- clarifying biological priorities
- normalizing user intent into config

Key modules:

- `module_preflight.R`
- `module_intake.R`
- `module_config.R`

### 2. Data and Validation Layer

Responsible for:

- reading count inputs
- reading metadata and annotation
- validating sample alignment and parameter policy

Key modules:

- `module_io.R`
- `module_validation.R`

### 3. Statistical Core

Responsible for:

- DESeq2 dataset construction
- fitting statistical models
- extracting differential results

Key modules:

- `module_design.R`
- `module_deseq2_engine.R`
- `module_contrasts.R`
- `module_deg_analysis.R`

### 4. Interpretation Layer

Responsible for:

- QC visualization
- trend clustering
- enrichment
- biological explanation

Key modules:

- `module_visualization.R`
- `module_tcseq_analysis.R`
- `module_enrichment.R`
- `module_biological_interpretation.R`
- `module_report_generation.R`

Key resources:

- `references/interpretation_rules.json`

### 5. Selection and Output Layer

Responsible for:

- variant comparison
- final selection
- manifest writing
- output organization

Key modules:

- `module_stability.R`
- `module_result_selection.R`
- `module_workflow_orchestrator.R`

## Planner/Agent Compatibility

The workflow should support:

- preflight-only mode
- validation-only mode
- DE-only mode
- full workflow mode
- report-only regeneration mode

To support that cleanly, module interfaces should remain functional and side-effect-minimal wherever possible.

## Interpretation Knowledge Layer

Biological explanation should not be hardcoded only in report templates. Use a separate interpretation knowledge layer so the skill can evolve without rewriting reporting code.

- `module_biological_interpretation.R` should resolve rule sets from study context, organism, and user-prioritized markers.
- `references/interpretation_rules.json` should hold reusable context tags, theme panels, and gene examples.
- Report generation should describe which rule panels were active for the current study so the interpretation remains auditable.
