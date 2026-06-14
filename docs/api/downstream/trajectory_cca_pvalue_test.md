# `cca_pvalue_test`

Permutation test for the CCA correlation. Shuffles the trajectory labels `num_simulations` times, recomputes CCA on the first two columns of the named DR, and builds an empirical null distribution of `|corr|` (the canonical correlation is sign-invariant). Returns the p-value (fraction of permutations whose correlation is at least the observed one) and writes a histogram of the null with the observed value marked.

**Source:** `sample_trajectory/CCA_test.py:519`

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
| `pseudo_adata` | AnnData | — | Sample-level AnnData whose `.obs` rows are samples and whose `.uns` holds the DR coordinates. |
| `column` | str | — | DR key in `.uns` (the single sample embedding, `"X_DR_sample"`); the first two columns are used. |
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
- `cca_pvalue_result_{column}.txt`

## Usage

```python
from sampledisco.sample_trajectory.CCA_test import cca_pvalue_test

pvalue = cca_pvalue_test(
    pseudo_adata=pseudo_adata,
    column="X_DR_sample",
    input_correlation=cca_score_a,
    output_directory="/results/rna",
    num_simulations=1000,
    trajectory_col="sev.level",
)
```
