# `cell_types_linux`

**Module**: `code/preparation/cell_type_gpu.py`

```python
cell_types_linux(
    anndata_cell,
    anndata_sample=None,
    cell_type_column="cell_type",
    existing_cell_types=False,
    n_target_clusters=None,
    umap=True,
    save=False,
    output_dir=None,
    defined_output_path=None,
    defined_sample_output_path=None,
    leiden_cluster_resolution=0.8,
    cell_embedding_column=None,
    cell_embedding_num_PCs=20,
    verbose=True,
    umap_plots=True,
    _recursion_depth=0,
)
```

Assigns or refines cell type labels with GPU-accelerated neighborhood/Leiden workflow.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `anndata_cell` | required | Cell-level AnnData. |
| `anndata_sample` | `None` | Optional sample-level AnnData to keep synchronized labels. |
| `cell_type_column` | `"cell_type"` | Input label column name if existing labels are reused. |
| `existing_cell_types` | `False` | Reuse existing labels instead of new Leiden clustering. |
| `n_target_clusters` | `None` | Target cluster count; enables iterative resolution adjustment and dendrogram aggregation. |
| `umap` | `True` | Compute UMAP embedding for visualization. |
| `save` | `False` | Save output AnnData files. |
| `output_dir` | `None` | Output directory when `save=True`. |
| `defined_output_path` | `None` | Explicit output path for cell-level object. |
| `defined_sample_output_path` | `None` | Explicit output path for sample-level object. |
| `leiden_cluster_resolution` | `0.8` | Leiden resolution seed value. |
| `cell_embedding_column` | `None` | Embedding key for clustering (`auto-detects` PCA vs LSI Harmony). |
| `cell_embedding_num_PCs` | `20` | PCs used for graph construction when PCA-based. |
| `verbose` | `True` | Print progress logs. |
| `umap_plots` | `True` | Save UMAP plot outputs when UMAP is computed. |
| `_recursion_depth` | `0` | Internal recursion counter for target-cluster mode. |

## Returns

- Updated `anndata_cell`
- Updated `anndata_sample` (if provided)

