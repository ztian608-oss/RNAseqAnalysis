# Parameter Policy

This skill is designed for professional RNA-seq differential expression analysis. Parameter choices must remain within common bioinformatics practice unless the user explicitly requests exploratory analysis and accepts that the result is non-standard.

## Allowed Thresholds

### log2 fold change

Allowed values:

- `1`
- `1.5`
- `2`

Interpretation:

- `1`: standard default, balances sensitivity and interpretability.
- `1.5`: moderate stringency.
- `2`: high stringency.

Forbidden examples:

- `0.3`
- `0.5`
- `0.58`
- `0.7`
- arbitrary decimals outside the allowed set

### Adjusted p-value

Allowed values:

- `0.05`
- `0.01`

Interpretation:

- `0.05`: default FDR threshold.
- `0.01`: stricter evidence threshold.

## TCseq Cluster Counts

Cluster-number comparison is allowed, but should remain interpretable. Recommended values:

- `4`
- `6`
- `8`

Do not generate large cluster grids by default. Use a small interpretable set and compare stability.

## Selection Principle

Final parameter selection should consider:

1. DEG stability
2. pathway stability
3. cluster stability when TCseq is used
4. alignment with the user’s biological question

Do not choose a final variant only because it returns the largest DEG count.
