# `cell_types`

Cell-type assignment via Leiden clustering on the sample-removed cell embedding (`Z_clust` by default). When `n_target_clusters` is set, the resolution is increased recursively until the Leiden output reaches the target; if it overshoots, a dendrogram-based aggregation merges the closest clusters (by centroid distance) until the count matches. Optionally computes and saves a UMAP. The resulting labels are written to `adata.obs["cell_type"]` as stringified integers.

**Source:** `preparation/cell_type_cpu.py:12` (GPU variant: `preparation/cell_type_gpu.py`, function `cell_types_gpu`)

## Signature

```python
def cell_types(
    anndata_cell,
    cell_type_column="cell_type",
    existing_cell_types=False,
    n_target_clusters=None,
    umap=True,
    save=False,
    output_dir=None,
    defined_output_path=None,
    leiden_cluster_resolution=0.8,
    cell_embedding_column=None,
    cell_embedding_num_PCs=20,
    verbose=True,
    umap_plots=True,
    _recursion_depth=0,
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `anndata_cell` | AnnData | — | Cell-level object from `preprocess` carrying `obsm["Z_clust"]` (sample-removed embedding). |
| `cell_type_column` | str | `"cell_type"` | Column name where labels are stored in `.obs`. |
| `existing_cell_types` | bool | `False` | If `True`, reuse labels already present in `cell_type_column` instead of re-clustering (and aggregate via dendrogram if above `n_target_clusters`). |
| `n_target_clusters` | int, optional | `None` | Target number of clusters. Triggers recursive resolution tuning + dendrogram aggregation. |
| `umap` | bool | `True` | Whether to compute a UMAP. |
| `save` | bool | `False` | Write the updated AnnData to disk. |
| `output_dir` | str, optional | `None` | Directory for the UMAP PNG, the `cell_type.csv`, and, with `save=True`, the updated h5ad. |
| `defined_output_path` | str, optional | `None` | Override the exact h5ad path that is written. |
| `leiden_cluster_resolution` | float | `0.8` | Starting Leiden resolution. |
| `cell_embedding_column` | str, optional | `None` | Which embedding to cluster on. When `None`, defaults to `"Z_clust"`. |
| `cell_embedding_num_PCs` | int | `20` | Number of components used in the neighborhood graph (RNA). |
| `verbose` | bool | `True` | Print progress. |
| `umap_plots` | bool | `True` | Save the UMAP visualization as a PNG. |
| `_recursion_depth` | int | `0` | Internal — do not set manually. |

## Returns

`AnnData` — the updated `anndata_cell` with `.obs["cell_type"]` populated and, if requested, `X_umap` in `.obsm`.

!!! note "API change"
    The previous two-argument form (`anndata_sample` plus `defined_sample_output_path`) and the `(anndata_cell, anndata_sample)` tuple return are gone. `cell_types` now operates on and returns a single cell-level AnnData; the sample-level object is materialized later by the sample-embedding step.

## Output files

- `{output_dir}/preprocess/umap_*.png` (when `umap_plots=True`, `umap=True`, and `output_dir` is set).
- `{output_dir}/preprocess/cell_type.csv` (cell-id → cell-type table, written whenever `output_dir` is set).
- `{output_dir}/preprocess/adata_preprocessed.h5ad` (or `defined_output_path`) when `save=True`.

## Usage

```python
from sampledisco.preparation.cell_type_cpu import cell_types

adata_cluster = cell_types(
    anndata_cell=adata_cluster,
    leiden_cluster_resolution=0.99,
    n_target_clusters=None,
    umap=True,
    save=True,
    output_dir="/results/rna",
)
```

!!! info "ATAC has its own function"
    ATAC cell typing is handled by `cell_types_atac` in `preparation/ATAC_cell_type.py` (GPU: `cell_types_atac_gpu`), which works on the ATAC diffusion-map / LSI embedding and uses a cosine neighborhood metric. The RNA `cell_types` above does detect ATAC inputs by the presence of `X_lsi` in `.obsm` and switches to a cosine metric for the neighborhood graph, but the dedicated `cell_types_atac` entry point is preferred for ATAC.
