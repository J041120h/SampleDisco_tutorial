# `find_optimal_cell_resolution_linux`

**Module**: `code/parameter_selection/gpu_optimal_resolution.py`

```python
find_optimal_cell_resolution_linux(
    adata_cell,
    adata_sample,
    output_dir,
    column,
    modality="rna",
    trajectory_col="sev.level",
    batch_col=None,
    sample_col="sample",
    celltype_col="cell_type",
    cell_embedding_column="X_pca_harmony",
    cell_embedding_num_pcs=20,
    n_hvg_features=2000,
    sample_embedding_dimension=10,
    harmony_for_proportion=True,
    preserve_cols_in_sample_embedding=None,
    n_cca_pcs=2,
    compute_corrected_pvalues=True,
    coarse_start=0.1,
    coarse_end=1.0,
    coarse_step=0.1,
    fine_range=0.02,
    fine_step=0.01,
    verbose=True,
)
```

GPU/Linux resolution search for RNA/ATAC cell type clustering.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata_cell` | required | Cell-level AnnData used for clustering. |
| `adata_sample` | required | Sample-level AnnData used for embedding/CCA scoring. |
| `output_dir` | required | Directory for optimization outputs. |
| `column` | required | DR key to optimize (`X_DR_expression` or `X_DR_proportion`). |
| `modality` | `"rna"` | Pipeline mode (`"rna"` or `"atac"`). |
| `trajectory_col` | `"sev.level"` | Trajectory supervision column for CCA scoring. |
| `batch_col` | `None` | Batch covariate(s) used in sample embedding. |
| `sample_col` | `"sample"` | Sample ID column. |
| `celltype_col` | `"cell_type"` | Cell type column name. |
| `cell_embedding_column` | `"X_pca_harmony"` | Cell embedding key used for Leiden clustering. |
| `cell_embedding_num_pcs` | `20` | Number of PCs for neighbor graph construction. |
| `n_hvg_features` | `2000` | HVG/HVF count used during sample embedding recomputation. |
| `sample_embedding_dimension` | `10` | Target sample DR dimension. |
| `harmony_for_proportion` | `True` | Harmony on proportion embedding branch. |
| `preserve_cols_in_sample_embedding` | `None` | Metadata columns preserved during embedding. |
| `n_cca_pcs` | `2` | CCA components used for trajectory score. |
| `compute_corrected_pvalues` | `True` | Compute corrected p-values for resolution scan. |
| `coarse_start` | `0.1` | Coarse search start resolution. |
| `coarse_end` | `1.0` | Coarse search end resolution. |
| `coarse_step` | `0.1` | Coarse search step size. |
| `fine_range` | `0.02` | Fine search range around coarse optimum. |
| `fine_step` | `0.01` | Fine search step size. |
| `verbose` | `True` | Print progress logs. |

## Returns

- `optimal_resolution`
- `results_df`

