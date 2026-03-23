# `replace_optimal_dimension_reduction`

**Module**: `code/parameter_selection/multi_omics_unify_optimal.py`

```python
replace_optimal_dimension_reduction(
    base_path,
    expression_resolution_dir=None,
    proportion_resolution_dir=None,
    pseudobulk_path=None,
    optimization_target="rna",
    verbose=True,
)
```

Replaces pseudobulk expression/proportion DR outputs with selected optimal-resolution results.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `base_path` | required | Base multi-omics output directory. |
| `expression_resolution_dir` | `None` | Custom path to expression optimization directory. |
| `proportion_resolution_dir` | `None` | Custom path to proportion optimization directory. |
| `pseudobulk_path` | `None` | Custom path to `pseudobulk_sample.h5ad`. |
| `optimization_target` | `"rna"` | Target modality used in optimization filenames (`rna` or `atac`). |
| `verbose` | `True` | Print replacement summary logs. |

## Returns

- Updated pseudobulk `AnnData` with replaced embeddings and expression matrix.

