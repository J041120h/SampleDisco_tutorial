# `cluster`

**Module**: `code/cluster.py`

```python
cluster(
    pseudobulk_adata,
    output_dir,
    number_of_clusters=5,
    use_expression=True,
    use_proportion=True,
    random_state=0,
)
```

K-means clustering for expression and proportion sample embeddings.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `pseudobulk_adata` | required | AnnData with `X_DR_expression`/`X_DR_proportion`. |
| `output_dir` | required | Output root (`sample_cluster/` created inside). |
| `number_of_clusters` | `5` | K in K-means. |
| `use_expression` | `True` | Cluster expression embedding. |
| `use_proportion` | `True` | Cluster proportion embedding. |
| `random_state` | `0` | Random seed for reproducibility. |

## Returns

- `expr_results`: `{sample_id: cluster}` or `None`
- `prop_results`: `{sample_id: cluster}` or `None`

