# `raisinfit`

**Module**: `code/sample_clustering/RAISIN.py`

```python
raisinfit(
    adata,
    sample_col,
    testtype="unpaired",
    group_col=None,
    individual_col=None,
    batch_col=None,
    sample_to_clade=None,
    custom_design=None,
    intercept=True,
    filtergene=False,
    filtergenequantile=0.5,
    n_jobs=None,
    verbose=True,
)
```

Fits RAISIN differential model for cluster/group comparisons.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Input AnnData for RAISIN fitting. |
| `sample_col` | required | Sample ID column. |
| `testtype` | `"unpaired"` | Design type (`paired`, `unpaired`, `continuous`, `custom`). |
| `group_col` | `None` | Group column in `obs` (preferred). |
| `individual_col` | `None` | Subject/individual ID column for paired designs. |
| `batch_col` | `None` | Batch column for optional ComBat correction. |
| `sample_to_clade` | `None` | Sample -> group mapping (used if `group_col` absent). |
| `custom_design` | `None` | Custom design matrices when `testtype="custom"`. |
| `intercept` | `True` | Include intercept in fixed-effects model. |
| `filtergene` | `False` | Enable low-expression gene filtering. |
| `filtergenequantile` | `0.5` | Quantile cutoff for gene filtering. |
| `n_jobs` | `None` | Parallel worker count (`None`/`-1` uses all cores). |
| `verbose` | `True` | Print progress logs. |

## Returns

- Fit dictionary consumed by `run_pairwise_tests`.

