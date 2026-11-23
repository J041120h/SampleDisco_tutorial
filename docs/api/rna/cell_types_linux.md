# `cell_types` (RNA and ATAC)

Cell-type assignment by Leiden clustering on the sample-removed cell embedding (`Z_comp` by default). The same function types RNA and ATAC cells; both wrappers call it.

- **Neighbor graph.** RNA: `sc.pp.neighbors` on the first `cell_embedding_num_PCs` columns of the embedding. ATAC (detected by `obsm['X_lsi']`): cosine metric on all columns of the embedding.
- **Clustering.** Leiden (igraph flavor, seed 42). Labels are written to `obs['cell_type']` as strings `"1"`, `"2"`, … .
- **Target count.** With `n_target_clusters`, the resolution is raised in steps of 0.5 (up to 5.0) until Leiden finds at least that many clusters; if it overshoots, cluster centroids are merged by average-linkage hierarchical clustering down to the target. The original labels are kept in `obs['cell_type_original']` and the merge map in `uns['cluster_mapping']`.
- **Existing labels.** With `existing_cell_types=True` and `cell_type_column` present in `.obs`, the labels are copied to `obs['cell_type']` instead of re-clustering. If there are more of them than `n_target_clusters`, they are merged down, with a `UserWarning` (the originals stay in `cell_type_column`; leave `n_target_clusters` at `None` to keep them all).

**Import:** `from sampledisco.preparation.cell_type_cpu import cell_types`

**Source:** `preparation/cell_type_cpu.py:13`. GPU variant: `cell_types_gpu` in `preparation/cell_type_gpu.py:14`, same signature.

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
| `anndata_cell` | AnnData | — | Cell-level object from [`preprocess`](preprocess_linux.md) (RNA) or [`preprocess` (ATAC)](../atac/preprocess_linux.md). Modified in place and returned. |
| `cell_type_column` | str | `"cell_type"` | Column **read** when `existing_cell_types=True`. Output labels always go to `obs['cell_type']`, whatever this is set to. |
| `existing_cell_types` | bool | `False` | Reuse the labels in `cell_type_column` instead of clustering. |
| `n_target_clusters` | int or None | `None` | Target number of cell types (see above). |
| `umap` | bool | `True` | Compute `obsm['X_umap']` (`min_dist=0.5`). |
| `save` | bool | `False` | Write the AnnData. Needs `output_dir` as well; see Output files. |
| `output_dir` | str or None | `None` | Root output directory. When set, `cell_type.csv` (and the UMAP plot) are written under `{output_dir}/preprocess/`. |
| `defined_output_path` | str or None | `None` | Exact h5ad path for `save=True`. Ignored unless `save=True` **and** `output_dir` is set. |
| `leiden_cluster_resolution` | float | `0.8` | Leiden resolution (the starting value when `n_target_clusters` is set). |
| `cell_embedding_column` | str or None | `None` | `.obsm` key to cluster on. `None` resolves `Z_comp`, then the legacy names `Z_clust`, `X_pca_harmony`, `X_lsi_harmony`, `X_glue_harmony`. |
| `cell_embedding_num_PCs` | int | `20` | Columns used for the RNA neighbor graph and the RNA merge step. Ignored for ATAC. |
| `verbose` | bool | `True` | Print progress. |
| `umap_plots` | bool | `True` | Save the UMAP plot (needs `umap=True` and `output_dir`). |
| `_recursion_depth` | int | `0` | Internal; do not set. |

## Returns

The same `AnnData` with `obs['cell_type']` (categorical), the neighbor graph in `.obsp`/`.uns['neighbors']`, and `obsm['X_umap']` when `umap=True`.

## Output files

With `output_dir` set:

- `{output_dir}/preprocess/cell_type.csv`: columns `cell_id`, `cell_type`.
- `{output_dir}/preprocess/umap_cell_type.png`: when `umap=True` and `umap_plots=True`.
- With `save=True`: the AnnData, gzip-compressed, at `defined_output_path` if given, otherwise `{output_dir}/preprocess/adata_preprocessed.h5ad` (overwriting the preprocess output, which is what the wrappers do).

## Usage

RNA, continuing from the [`preprocess`](preprocess_linux.md) example:

```python
import anndata as ad
from sampledisco.preparation.cell_type_cpu import cell_types

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
adata = cell_types(
    adata,
    n_target_clusters=8,
    umap=True,
    save=True,
    output_dir="api_demo/rna",
)
print(adata.obs["cell_type"].value_counts().sort_index())
```

**Writes** → `api_demo/rna/preprocess/cell_type.csv`, `umap_cell_type.png`, and the updated `adata_preprocessed.h5ad`. Leiden at resolution 0.8 finds 18 clusters, which are merged to 8:

```text
[cell_types] Found 18 clusters
[cell_types] Over-shot target; recursing with dendrogram aggregation
  [cell_types] Aggregating 18 → 8 using dendrogram
...
cell_type
1      29
2    1075
3    7949
4    7415
5    9925
6    1974
7     216
8    1097
```

Merging by centroid distance can leave very uneven types, as here. The cluster count can differ by one or two between machines (an earlier run found 19). About 75 s (most of it UMAP) and 2.2 GB.

ATAC, continuing from the [ATAC `preprocess`](../atac/preprocess_linux.md) example (the same call; the cosine graph is chosen automatically):

```python
import anndata as ad
from sampledisco.preparation.cell_type_cpu import cell_types

adata_atac = ad.read_h5ad("api_demo/atac/preprocess/adata_preprocessed.h5ad")
adata_atac = cell_types(
    adata_atac,
    n_target_clusters=8,
    umap=False,
    save=True,
    output_dir="api_demo/atac",
)
print(adata_atac.obs["cell_type"].value_counts().sort_index())
```

**Writes** → `api_demo/atac/preprocess/cell_type.csv` and the updated `adata_preprocessed.h5ad`. Leiden finds exactly 8 clusters here, so nothing is merged (an earlier run on another node found 9 and merged them to 8). About 40 s, 3.2 GB:

```text
cell_type
1    4902
2    7150
3    2280
4     999
5    4832
6    3317
7    1419
8     294
```

!!! warning "`cell_types_atac` is not used by the pipeline"
    `sampledisco.preparation.ATAC_cell_type.cell_types_atac` (GPU: `cell_types_atac_gpu`) is an older ATAC-only routine with a different signature (`Save=`, `use_rep=`, `cluster_resolution=`). The ATAC wrapper calls the `cell_types` documented here. Two quirks of `cell_types_atac`: with `Save=True` it writes `adata_sample.h5ad` to `{output_dir}/preprocess/` when it also draws its UMAP plots (`generate_plots=True` and a UMAP exists), and to `{output_dir}/` otherwise; and it uses a cosine neighbor graph only for existing labels, not for fresh clustering.

!!! note "Wrapper behaviour"
    The wrappers call `cell_types(..., cell_type_column=<*_celltype_col>, save=True, output_dir=<modality output dir>)`. They then switch to `obs['cell_type']` for the sample embedding and the downstream steps (and return the column used as `results['rna_results']['celltype_col']`), so the new or merged labels are used whatever `rna_celltype_col` / `atac_celltype_col` is. GPU (`cell_types_gpu`) is used only when the config sets `use_gpu: true` and RAPIDS imports.
