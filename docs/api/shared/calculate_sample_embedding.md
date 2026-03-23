# `calculate_sample_embedding`

**Module**: `code/sample_embedding/calculate_sample_embedding.py`

```python
calculate_sample_embedding(
    adata,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="./",
    sample_hvg_number=2000,
    combat_timeout=None,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    preserve_cols_in_sample_embedding=None,
    use_gpu=False,
    atac=False,
    save=True,
    verbose=True,
)
```

Builds sample-level embeddings from single-cell data.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Input single-cell AnnData. |
| `sample_col` | `"sample"` | Sample identifier column in `adata.obs`. |
| `celltype_col` | `"cell_type"` | Cell type label column in `adata.obs`. |
| `batch_col` | `None` | Batch column(s) for correction. |
| `output_dir` | `"./"` | Output directory for pseudobulk/embedding files. |
| `sample_hvg_number` | `2000` | Number of variable features per cell type for pseudobulk. |
| `combat_timeout` | `None` | Timeout for CPU ComBat stage. |
| `n_expression_components` | `10` | Reduced dimensions for expression embedding. |
| `n_proportion_components` | `10` | Reduced dimensions for proportion embedding. |
| `harmony_for_proportion` | `True` | Apply Harmony to proportion embedding. |
| `preserve_cols_in_sample_embedding` | `None` | Metadata columns preserved through embedding pipeline. |
| `use_gpu` | `False` | Use GPU pseudobulk implementation. |
| `atac` | `False` | Switches behavior for ATAC-style representation (`True` for ATAC). |
| `save` | `True` | Save generated artifacts to disk. |
| `verbose` | `True` | Print progress logs. |

## Returns

- `pseudobulk_result_dict`
- `sample_embedding_adata` (AnnData with sample embeddings)

