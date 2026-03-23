# `CCA_Call`

**Module**: `code/sample_trajectory/CCA.py`

```python
CCA_Call(
    adata,
    output_dir=None,
    trajectory_col="sev.level",
    n_components=2,
    auto_select_best_2pc=True,
    verbose=False,
    show_sample_labels=False,
)
```

Runs supervised CCA trajectory analysis on expression and proportion DR spaces.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Sample-level AnnData with DR coordinates. |
| `output_dir` | `None` | Output root (`CCA/` is created under this path). |
| `trajectory_col` | `"sev.level"` | Supervision column in `adata.obs`. |
| `n_components` | `2` | Number of components used for CCA/PCA views. |
| `auto_select_best_2pc` | `True` | Auto-select best 2 PC combination for plotting. |
| `verbose` | `False` | Print progress logs. |
| `show_sample_labels` | `False` | Draw sample names on plots. |

## Returns

- `(proportion_score, expression_score, proportion_pseudotime, expression_pseudotime)`

