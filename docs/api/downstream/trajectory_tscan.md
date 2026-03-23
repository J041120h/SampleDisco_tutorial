# `TSCAN`

**Module**: `code/sample_trajectory/TSCAN.py`

```python
TSCAN(
    AnnData_sample,
    column,
    n_clusters=None,
    output_dir="./",
    grouping_columns=None,
    verbose=False,
    origin=None,
    pseudotime_mode="rank",
)
```

Unsupervised sample trajectory inference using TSCAN-like GMM + MST pipeline.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `AnnData_sample` | required | Sample-level AnnData with DR coordinates. |
| `column` | required | DR key used as trajectory input. |
| `n_clusters` | `None` | Number of GMM clusters; auto-BIC when `None`. |
| `output_dir` | `"./"` | Output root (`TSCAN/` created inside). |
| `grouping_columns` | `None` | Metadata columns used in grouping visualizations. |
| `verbose` | `False` | Print detailed progress logs. |
| `origin` | `None` | Optional starting cluster index for pseudotime path. |
| `pseudotime_mode` | `"rank"` | Pseudotime scaling mode (`rank` or `distance`). |

## Returns

- Dict containing clusters, MST path, projections, and pseudotime outputs.

