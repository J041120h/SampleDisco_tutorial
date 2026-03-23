# `cca_pvalue_test`

**Module**: `code/sample_trajectory/CCA_test.py`

```python
cca_pvalue_test(
    pseudo_adata,
    column,
    input_correlation,
    output_directory,
    num_simulations=1000,
    trajectory_col="sev.level",
    verbose=True,
)
```

Permutation-based p-value test for observed CCA trajectory correlation.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `pseudo_adata` | required | Sample-level AnnData containing DR coordinates. |
| `column` | required | Coordinate key (e.g., `X_DR_expression`). |
| `input_correlation` | required | Observed CCA correlation to test. |
| `output_directory` | required | Output root (`CCA_test/` created inside). |
| `num_simulations` | `1000` | Number of permutation simulations. |
| `trajectory_col` | `"sev.level"` | Trajectory label column in `obs`. |
| `verbose` | `True` | Print timing/progress logs. |

## Returns

- `p_value` (float)

