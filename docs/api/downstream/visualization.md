# `visualization`

General diagnostic plots, each behind a flag:

1. **Cell-type dendrogram** (`plot_dendrogram_flag`): the mean cell embedding of each cell type (`Z_comp`, else a legacy composition key, else `X_glue`), Euclidean distances, Ward linkage. The cell-type column is the first of `cell_type`, `celltype`, `cluster`, `leiden`, `seurat_clusters` found in `AnnData_cell.obs`.
2. **"Proportion" scatter** (`plot_cell_type_proportions_pca_flag`): samples on PC1 vs PC2 of `X_DR_sample`, colored by each grouping column (`RdYlBu` colormap).
3. **"Expression" scatter** (`plot_cell_type_expression_umap_flag`): the **same** PC1 vs PC2 scatter of `X_DR_sample`, `viridis` colormap. Despite the flag names, neither plot is a PCA of proportions or a UMAP of expression; both names survive from the old two-embedding layout.

**Import:** `from sampledisco.visualization.visualization_other import visualization`

**Source:** `visualization/visualization_other.py:71`

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
| `AnnData_cell` | AnnData | — | Cell-level AnnData (cell types and the cell embedding, for the dendrogram). |
| `pseudobulk_anndata` | AnnData | — | Sample-level AnnData ([`build_sample_adata`](../shared/build_sample_adata.md)) with `X_DR_sample` and the grouping columns in `.obs`. |
| `output_dir` | str | — | Output folder; see Output files for the `visualization/` subfolder rule. |
| `grouping_columns` | list or None | `None` | `.obs` columns of `pseudobulk_anndata`, one scatter per column. Required (non-empty) for the two scatter flags. Missing columns raise `ValueError`. |
| `age_bin_size` | int or None | `None` | Only adds `age_column` to the columns checked to exist; no binning is done. |
| `age_column` | str | `"age"` | See `age_bin_size`. |
| `verbose` | bool | `True` | Print progress. |
| `plot_dendrogram_flag` | bool | `True` | Write the dendrogram. |
| `plot_cell_type_proportions_pca_flag` | bool | `False` | Write the `proportion_embedding_*` scatters. |
| `plot_cell_type_expression_umap_flag` | bool | `False` | Write the `expression_embedding_*` scatters. |

## Returns

`None`.

## Output files

With `grouping_columns` set, the files go to `{output_dir}/visualization/`; with `grouping_columns=None`, directly to `{output_dir}/`.

- `cell_type_dendrogram.png`
- `proportion_embedding_{column}.png`, one per grouping column
- `expression_embedding_{column}.png`, one per grouping column

## Usage

Continuing from the RNA chain:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.visualization.visualization_other import visualization

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

visualization(
    adata,
    pseudo,
    "api_demo/rna",
    grouping_columns=["sev.level"],
    plot_dendrogram_flag=True,
    plot_cell_type_proportions_pca_flag=True,
)
```

**Writes** → `api_demo/rna/visualization/cell_type_dendrogram.png` (8 cell types) and `api_demo/rna/visualization/proportion_embedding_sev.level.png`. About 2 s.
