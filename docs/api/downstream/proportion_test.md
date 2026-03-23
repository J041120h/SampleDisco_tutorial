# `proportion_test`

**Module**: `code/sample_clustering/proportion_test.py`

```python
proportion_test(
    adata,
    sample_col,
    group_col=None,
    sample_to_clade=None,
    celltype_col="celltype",
    output_dir=None,
    verbose=True,
)
```

Tests cell type composition differences across sample groups.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Cell-level AnnData. |
| `sample_col` | required | Sample ID column. |
| `group_col` | `None` | Group column in `obs` (preferred when available). |
| `sample_to_clade` | `None` | Sample -> group mapping used when `group_col` is absent. |
| `celltype_col` | `"celltype"` | Cell type column in `obs`. |
| `output_dir` | `None` | Output path for plots/tables. |
| `verbose` | `True` | API compatibility flag (minimal runtime effect). |

## Notes

- Provide either `group_col` or `sample_to_clade`.
- If both are provided, `group_col` is used.

