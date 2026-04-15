# `preprocess_linux` (ATAC)

End-to-end preprocessing for scATAC-seq on the GPU/Linux code path. The function reads a peak-count `.h5ad`, merges metadata, applies quality-control filters, runs TF-IDF normalization, projects via LSI, and applies Harmony for batch correction on the clustering track. It also optionally runs scDblFinder-style doublet detection and can drop the first LSI component (which typically reflects sequencing depth). Returns the same two-AnnData pair as the RNA version: a Harmony-integrated object for clustering and a minimally processed object for sample-level analysis.

**Source:** `preparation/atac_preprocess_gpu.py:204`

## Signature

```python
def preprocess_linux(
    h5ad_path,
    sample_meta_path,
    output_dir,
    cell_meta_path=None,
    sample_column="sample",
    sample_level_batch_key="batch",
    cell_embedding_num_PCs=50,
    num_harmony_iterations=30,
    num_cell_hvfs=50000,
    min_cells=1,
    min_features=2000,
    max_features=15000,
    min_cells_per_sample=1,
    exclude_features=None,
    cell_level_batch_key=None,
    doublet_detection=True,
    tfidf_scale_factor=1e4,
    log_transform=True,
    drop_first_lsi=True,
    verbose=True,
) -> Tuple[AnnData, AnnData]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `h5ad_path` | str | — | Cell × peak `.h5ad`. |
| `sample_meta_path` | str | — | Per-sample metadata CSV. |
| `output_dir` | str | — | Writes to `{output_dir}/preprocess/`. |
| `cell_meta_path` | str, optional | `None` | Per-cell metadata CSV merged into `.obs`. |
| `sample_column` | str | `"sample"` | Sample identifier column. |
| `sample_level_batch_key` | str or list, optional | `"batch"` | Sample-level Harmony batch. Set `None` to skip. |
| `cell_embedding_num_PCs` | int | `50` | Number of LSI components. |
| `num_harmony_iterations` | int | `30` | Max Harmony iterations. |
| `num_cell_hvfs` | int | `50000` | Number of highly variable features for the clustering track. |
| `min_cells` | int | `1` | Drop features seen in fewer than this many cells. |
| `min_features` | int | `2000` | Drop cells with fewer than this many features. |
| `max_features` | int | `15000` | Drop cells with more than this many features. |
| `min_cells_per_sample` | int | `1` | Drop samples with fewer than this many surviving cells. |
| `exclude_features` | list, optional | `None` | Feature IDs to drop before HVF selection. |
| `cell_level_batch_key` | list, optional | `None` | Cell-level Harmony batch keys; sample id is always included. |
| `doublet_detection` | bool | `True` | Run doublet detection before clustering. |
| `tfidf_scale_factor` | float | `1e4` | Scaling factor inside TF-IDF normalization. |
| `log_transform` | bool | `True` | Apply `log1p` after TF-IDF. |
| `drop_first_lsi` | bool | `True` | Discard LSI component 1 (typically depth-correlated). |
| `verbose` | bool | `True` | Print progress. |

## Returns

`Tuple[AnnData, AnnData]` — `(adata_cluster, adata_sample)`:

- `adata_cluster` carries `X_lsi_harmony` in `.obsm`.
- `adata_sample` preserves raw peak counts for downstream pseudobulk.

## Output files

- `{output_dir}/preprocess/adata_cell.h5ad`
- `{output_dir}/preprocess/adata_sample.h5ad`

## Usage

```python
from genodistance.preparation import preprocess_linux  # ATAC build

adata_cluster, adata_sample = preprocess_linux(
    h5ad_path="/data/test_ATAC.h5ad",
    sample_meta_path="/data/sample_meta.csv",
    output_dir="/results/atac",
    cell_embedding_num_PCs=50,
    num_cell_hvfs=50000,
    doublet_detection=True,
    drop_first_lsi=True,
)
```

!!! note "Cell typing is shared with RNA"
    ATAC does not have its own `cell_types_linux` — the RNA version auto-detects `X_lsi_harmony` and switches to cosine neighborhood metrics. See [`cell_types_linux`](../rna/cell_types_linux.md).
