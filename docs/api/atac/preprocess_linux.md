# `preprocess` (ATAC)

End-to-end preprocessing for scATAC-seq. The function reads a peak-count `.h5ad`, merges metadata, applies quality-control filters, runs TF-IDF normalization, projects via LSI, and applies a two-pass Harmony correction. It also optionally runs Scrublet doublet detection and can drop the first LSI component (which typically reflects sequencing depth). Like the RNA version, it writes a single preprocessed AnnData carrying both a sample-removed clustering embedding (`obsm['Z_clust']`) and a sample-preserved displacement embedding (`obsm['Z_rmd']`), and returns that AnnData.

The CPU implementation lives in `preparation/atac_preprocess_cpu.py`; the GPU-accelerated variant `preprocess_gpu` lives in `preparation/atac_preprocess_gpu.py` and activates only when the RAPIDS stack is importable.

## Signature

```python
def preprocess(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
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
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `h5ad_path` | str | — | Cell × peak `.h5ad`. |
| `sample_meta_path` | str | — | Per-sample metadata CSV. |
| `output_dir` | str | — | Writes to `{output_dir}/preprocess/`. |
| `sample_column` | str | `"sample"` | Sample identifier column. |
| `cell_meta_path` | str, optional | `None` | Per-cell metadata CSV (indexed by `barcode`) merged into `.obs`. |
| `sample_level_batch_key` | str or list, optional | `None` | Sample-level Harmony batch. Leave `None` to skip. |
| `cell_embedding_num_PCs` | int | `50` | Number of LSI components. |
| `num_harmony_iterations` | int | `30` | Max Harmony iterations. |
| `num_cell_hvfs` | int | `50000` | Number of highly variable features for the clustering track. |
| `min_cells` | int | `1` | Drop features seen in fewer than this many cells. |
| `min_features` | int | `2000` | Drop cells with fewer than this many features. |
| `max_features` | int | `15000` | Drop cells with more than this many features. |
| `min_cells_per_sample` | int | `1` | Drop samples with fewer than this many surviving cells. |
| `exclude_features` | list, optional | `None` | Feature IDs to drop before HVF selection. |
| `cell_level_batch_key` | list, optional | `None` | Cell-level Harmony batch keys; sample id is always included. |
| `doublet_detection` | bool | `True` | Run Scrublet doublet detection before clustering. |
| `tfidf_scale_factor` | float | `1e4` | Scaling factor inside TF-IDF normalization. |
| `log_transform` | bool | `True` | Apply `log1p` after TF-IDF. |
| `drop_first_lsi` | bool | `True` | Discard LSI component 1 (typically depth-correlated). |
| `verbose` | bool | `True` | Print progress. |

## Returns

`AnnData` — a single preprocessed object. The raw peak counts are kept in `.layers['counts']`, `.X` holds the TF-IDF + `log1p` normalized matrix, and `.obsm` carries:

- `X_lsi` — LSI on the HVF subset (with the first component dropped when `drop_first_lsi=True`).
- `Z_clust` — sample-removed Harmony embedding (clustering / composition view).
- `Z_rmd` — sample-preserved Harmony embedding (RMD / displacement view).

## Output files

- `{output_dir}/preprocess/adata_preprocessed.h5ad`

## Usage

```python
from sampledisco.preparation.atac_preprocess_cpu import preprocess  # GPU: atac_preprocess_gpu.preprocess_gpu

adata = preprocess(
    h5ad_path="data/test_ATAC.h5ad",
    sample_meta_path=None,
    output_dir="sampledisco_demo_output/atac",
    cell_embedding_num_PCs=50,
    num_cell_hvfs=50000,
    doublet_detection=True,
    drop_first_lsi=True,
)
```

!!! note "Prefer the config-driven wrapper"
    These low-level functions are the building blocks. The supported entry point is the YAML-config wrapper, which runs ATAC preprocessing → cell typing → sample embedding → downstream analysis in one call: `sampledisco -m complex --config <config.yaml>` (set `run_atac_pipeline: true`).

!!! note "ATAC cell typing"
    ATAC has its own cell-typing routine, `cell_types_atac` (`from sampledisco.preparation.ATAC_cell_type import cell_types_atac`), which builds a dendrogram and differential peaks on the ATAC DR embedding. The shared RNA typing routine `cell_types` operates on `Z_clust`. See [`cell_types_linux`](../rna/cell_types_linux.md).
