# `cca_pvalue_test`

Permutation test for the CCA correlation. Shuffles the trajectory labels `num_simulations` times, recomputes CCA on the first two components of the named DR, and builds an empirical null distribution. Returns the p-value (fraction of permutations whose correlation exceeds the observed one) and writes a density plot with the observed value marked.

**Source:** `sample_trajectory/CCA_test.py:599`

## Signature

```python
def cca_pvalue_test(
    pseudo_adata: AnnData,
    column: str,
    input_correlation: float,
    output_directory: str,
    num_simulations: int = 1000,
    trajectory_col: str = "sev.level",
    verbose: bool = True,
) -> float
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudo_adata` | AnnData | — | Sample-level pseudobulk AnnData with DR results in `.uns`. |
| `column` | str | — | DR key in `.uns` (e.g. `"X_DR_expression"` or `"X_DR_proportion"`). |
| `input_correlation` | float | — | Observed CCA correlation to test against (from [`CCA_Call`](trajectory_cca_call.md)). |
| `output_directory` | str | — | Writes to `{output_directory}/CCA_test/`. |
| `num_simulations` | int | `1000` | Number of label permutations for the null distribution. |
| `trajectory_col` | str | `"sev.level"` | Phenotype column in `pseudo_adata.obs`. |
| `verbose` | bool | `True` | Print timing. |

## Returns

`float` — empirical p-value.

## Output files

Under `{output_directory}/CCA_test/`:

- `cca_pvalue_distribution_{column}.png`
- `cca_pvalue_results_{column}.csv`

## Usage

```python
from genodistance.sample_trajectory import cca_pvalue_test

pvalue = cca_pvalue_test(
    pseudo_adata=pseudo_adata,
    column="X_DR_expression",
    input_correlation=expr_score,
    output_directory="/results/rna",
    num_simulations=1000,
    trajectory_col="sev.level",
)
```
