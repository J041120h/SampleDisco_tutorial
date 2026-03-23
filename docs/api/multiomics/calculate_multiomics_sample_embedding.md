# `calculate_multiomics_sample_embedding`

**Module**: `code/sample_embedding/calculate_multiomics_sample_embedding.py`

```python
calculate_multiomics_sample_embedding(
    adata,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="./",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    preserve_cols_in_sample_embedding=None,
    use_gpu=False,
    atac=False,
    save=True,
    verbose=True,
    hvg_modality="RNA",
    modality_col="modality",
)
```

Computes pseudobulk and sample embeddings for integrated multi-omics data.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Integrated AnnData input. |
| `sample_col` | `"sample"` | Sample ID column. |
| `celltype_col` | `"cell_type"` | Cell type label column. |
| `batch_col` | `None` | Batch covariate(s) for correction. |
| `output_dir` | `"./"` | Output root for pseudobulk and embedding files. |
| `sample_hvg_number` | `2000` | Number of features selected for pseudobulk. |
| `n_expression_components` | `10` | Expression embedding dimension. |
| `n_proportion_components` | `10` | Proportion embedding dimension. |
| `harmony_for_proportion` | `True` | Apply Harmony on proportion embedding. |
| `preserve_cols_in_sample_embedding` | `None` | Metadata columns preserved in processing. |
| `use_gpu` | `False` | Use GPU pseudobulk path. |
| `atac` | `False` | ATAC mode switch in shared routines. |
| `save` | `True` | Save generated files. |
| `verbose` | `True` | Print progress logs. |
| `hvg_modality` | `"RNA"` | Modality used for HVG selection. |
| `modality_col` | `"modality"` | Modality column in `obs`. |

## Returns

- `pseudobulk_result_dict`
- `sample_embedding_adata`

