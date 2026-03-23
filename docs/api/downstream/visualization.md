# `visualization`

**Module**: `code/visualization/visualization_other.py`

```python
visualization(
    AnnData_cell,
    pseudobulk_anndata,
    output_dir,
    grouping_columns=None,
    age_bin_size=None,
    age_column="age",
    verbose=True,
    plot_dendrogram_flag=True,
    plot_cell_type_proportions_pca_flag=False,
    plot_cell_type_expression_umap_flag=False,
)
```

Runs downstream visualization submodules conditionally.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `AnnData_cell` | required | Cell-level AnnData (for dendrogram). |
| `pseudobulk_anndata` | required | Sample-level AnnData for embedding plots. |
| `output_dir` | required | Output root for visualization files. |
| `grouping_columns` | `None` | Metadata columns used for grouped plots. |
| `age_bin_size` | `None` | Optional age binning size in preprocessing step. |
| `age_column` | `"age"` | Age column name for age-derived grouping. |
| `verbose` | `True` | Print progress logs. |
| `plot_dendrogram_flag` | `True` | Create cell type dendrogram. |
| `plot_cell_type_proportions_pca_flag` | `False` | Plot proportion embedding by grouping columns. |
| `plot_cell_type_expression_umap_flag` | `False` | Plot expression embedding by grouping columns. |

