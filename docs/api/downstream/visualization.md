# `visualization`

Bundle of general-purpose diagnostic plots. Given the cell-level AnnData and the sample-level pseudobulk, it can emit three categories of figure (each individually toggleable):

1. **Cell-type dendrogram** — hierarchical clustering of cell types based on their pseudobulk expression profiles.
2. **Cell-type proportion PCA** — PCA of per-sample cell-type proportions with samples colored by a user-specified grouping column.
3. **Cell-type expression UMAP** — UMAP of samples on the pseudobulk expression embedding.

Supply one or more `grouping_columns` (and optionally `age_bin_size` for age-based faceting) to control how samples are colored.

**Source:** `visualization/visualization_other.py:73`

## Signature

```python
def visualization(
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

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `AnnData_cell` | AnnData | — | Cell-level AnnData carrying the cell-type labels. |
| `pseudobulk_anndata` | AnnData | — | Sample-level pseudobulk with DR in `.obsm`. |
| `output_dir` | str | — | Writes plots under a grouping-aware subdirectory of `{output_dir}/visualization/`. |
| `grouping_columns` | list, optional | `None` | `.obs` columns used for sample coloring. Required for proportion PCA / expression UMAP plots. |
| `age_bin_size` | int, optional | `None` | Optional age bin width; enables age-as-category coloring. |
| `age_column` | str | `"age"` | Column holding sample age. |
| `verbose` | bool | `True` | Print progress. |
| `plot_dendrogram_flag` | bool | `True` | Emit the cell-type dendrogram PNG. |
| `plot_cell_type_proportions_pca_flag` | bool | `False` | Emit proportion PCAs per grouping column. |
| `plot_cell_type_expression_umap_flag` | bool | `False` | Emit expression UMAPs per grouping column. |

## Returns

`None`. Side effects only.

## Output files

Under `{output_dir}/visualization/`:

- `cell_type_dendrogram.png`
- `proportion_pca_{grouping_column}.png` (per flag)
- `expression_umap_{grouping_column}.png` (per flag)

## Usage

```python
from genodistance.visualization import visualization

visualization(
    AnnData_cell=adata_cell,
    pseudobulk_anndata=pseudo_adata,
    output_dir="/results/rna",
    grouping_columns=["sev.level"],
    plot_dendrogram_flag=True,
    plot_cell_type_proportions_pca_flag=True,
    plot_cell_type_expression_umap_flag=False,
)
```
