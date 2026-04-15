# `cell_types_linux`

Cell-type assignment via Leiden clustering on the Harmony-integrated cell embedding (`X_pca_harmony` for RNA or `X_lsi_harmony` for ATAC — auto-detected from `adata.obsm`). When `n_target_clusters` is set, the resolution is adjusted recursively until the Leiden output matches the target; if it overshoots, a dendrogram-based aggregation merges the smallest clusters until the count matches. Optionally computes and saves a UMAP. The resulting labels are written to `adata.obs["cell_type"]` as stringified integers.

**Source:** `preparation/cell_type_gpu.py:13`

## Signature

```python
def cell_types_linux(
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
) -> Tuple[AnnData, AnnData]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `anndata_cell` | AnnData | — | Cell-level object from `preprocess_linux` with a Harmony-corrected embedding in `.obsm`. |
| `anndata_sample` | AnnData, optional | `None` | Sample-level object; labels are propagated to it when provided. |
| `cell_type_column` | str | `"cell_type"` | Column name where labels are stored in `.obs`. |
| `existing_cell_types` | bool | `False` | If `True`, reuse labels already present in `cell_type_column` instead of re-clustering. |
| `n_target_clusters` | int, optional | `None` | Target number of clusters. Triggers recursive resolution tuning + dendrogram aggregation. |
| `umap` | bool | `True` | Whether to compute a UMAP. |
| `save` | bool | `False` | Write updated AnnData objects to disk. |
| `output_dir` | str, optional | `None` | Directory for the UMAP PNG and, with `save=True`, the updated h5ads. |
| `defined_output_path` / `defined_sample_output_path` | str, optional | `None` | Override the exact h5ad paths that are written. |
| `leiden_cluster_resolution` | float | `0.8` | Starting Leiden resolution. |
| `cell_embedding_column` | str, optional | `None` | Which embedding to cluster on. When `None`, prefers `X_lsi_harmony` if present (ATAC), else `X_pca_harmony` (RNA). |
| `cell_embedding_num_PCs` | int | `20` | Number of components used in the neighborhood graph. |
| `verbose` | bool | `True` | Print progress. |
| `umap_plots` | bool | `True` | Save the UMAP visualization as a PNG. |
| `_recursion_depth` | int | `0` | Internal — do not set manually. |

## Returns

`Tuple[AnnData, AnnData]` — updated `(anndata_cell, anndata_sample)` with `.obs["cell_type"]` populated and, if requested, `X_umap` in `.obsm`.

## Output files

- `{output_dir}/preprocess/umap_*.png` (when `umap_plots=True`).
- Updated h5ad files when `save=True`.

## Usage

```python
from genodistance.preparation import cell_types_linux

adata_cluster, adata_sample = cell_types_linux(
    anndata_cell=adata_cluster,
    anndata_sample=adata_sample,
    leiden_cluster_resolution=0.99,
    n_target_clusters=None,
    umap=True,
    save=True,
    output_dir="/results/rna",
)
```

!!! info "ATAC reuses this function"
    The ATAC pipeline imports the same `cell_types_linux`. It detects `X_lsi_harmony` in `.obsm` and automatically switches to the ATAC-appropriate neighborhood metric.
