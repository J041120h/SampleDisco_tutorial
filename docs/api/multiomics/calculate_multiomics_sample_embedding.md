# `compute_sample_embedding` (multi-omics)

!!! warning "Renamed and unified"
    There is no longer a separate `calculate_multiomics_sample_embedding` function. The unified [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) handles RNA, ATAC, **and** multi-omics from a single entry point. For multi-omics, pass `modality_col="modality"`; the function automatically tags units as `<sample>_RNA` / `<sample>_ATAC` and builds one shared sample-level embedding across modalities. The OLD two-key return (`X_DR_expression`, `X_DR_proportion`) is gone — the current pipeline produces a **single** key, `uns['X_DR_sample']`.

Lifts the integrated RNA + ATAC cell-level AnnData into a sample-level embedding. It combines multi-resolution cell-type **composition** blocks (computed on the sample-removed `Z_clust` view) with an **RMD displacement** block (computed on the sample-preserved `Z_rmd` view), then inverse-variance weights, Frobenius-stacks, PCA-reduces, and Harmony-corrects them at the sample level. The result is written in place to `adata.uns['X_DR_sample']` and the modified `AnnData` is returned.

**Source:** `sample_embedding/__init__.py:17` (CPU/GPU dispatcher) → `sample_embedding/sample_embedding.py:130`

## Signature

```python
def compute_sample_embedding(
    adata: AnnData,
    output_dir: str,
    *,
    use_gpu: bool = False,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    cluster_emb_key: str = "Z_clust",
    rmd_emb_key: Optional[str] = None,
    modality_col: Optional[str] = None,
    batch_col: Optional[Union[str, List[str]]] = None,
    medium_K: int = 120,
    fine_K: int = 300,
    rmd_dim_per_cluster: int = 8,
    use_clr: bool = False,
    use_rmd: bool = True,
    block_weights: Optional[List[float]] = None,
    rmd_weight: float = 0.60,
    pca_components: int = 10,
    batch_method: str = "harmony",
    save: bool = True,
    verbose: bool = True,
    seed: int = 42,
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Integrated cell-level AnnData with `cell_type` labels, `obsm['Z_clust']` (sample-removed) and `obsm['Z_rmd']` (sample-preserved, primary `X_glue` for multi-omics), and a `modality_col`. |
| `output_dir` | str | — | Directory for saved outputs. |
| `use_gpu` | bool | `False` | Dispatch to the RAPIDS-accelerated GPU implementation when the RAPIDS stack is importable; otherwise the CPU path is used. |
| `sample_col` | str | `"sample"` | Sample identifier column. |
| `celltype_col` | str | `"cell_type"` | Cell-type column. |
| `cluster_emb_key` | str | `"Z_clust"` | obsm key for the sample-removed (clustering / composition) embedding. |
| `rmd_emb_key` | str, optional | `None` | obsm key for the sample-preserved (RMD displacement) embedding; resolved automatically (e.g. to `Z_rmd`, legacy `Z_cmd`) when `None`. |
| `modality_col` | str, optional | `None` | Column identifying modality per cell. **Set to `"modality"` for multi-omics** so units are split into `<sample>_RNA` / `<sample>_ATAC`. |
| `batch_col` | str or list, optional | `None` | Sample-level batch covariate(s) for the sample-level Harmony correction. |
| `medium_K` | int | `120` | Number of clusters for the medium-resolution composition block. |
| `fine_K` | int | `300` | Number of clusters for the fine-resolution composition block. |
| `rmd_dim_per_cluster` | int | `8` | RMD displacement dimensions retained per coarse cell type. |
| `use_clr` | bool | `False` | CLR-transform the composition blocks. |
| `use_rmd` | bool | `True` | Include the RMD displacement block. |
| `block_weights` | list, optional | `None` | Explicit per-block weights; inverse-variance weighting is used when `None`. |
| `rmd_weight` | float | `0.60` | Relative weight of the RMD block versus the composition blocks. |
| `pca_components` | int | `10` | Number of PCs in the final sample embedding. |
| `batch_method` | str | `"harmony"` | Sample-level batch-correction method. |
| `save` | bool | `True` | Write outputs to `output_dir`. |
| `verbose` | bool | `True` | Print progress. |
| `seed` | int | `42` | Random seed. |

## Returns

`AnnData` — the **same** cell-level AnnData, modified in place, now carrying:

- `adata.uns['X_DR_sample']` — the final sample embedding as a pandas `DataFrame` (units × PCs). For multi-omics, units are `<sample>_RNA` / `<sample>_ATAC`.
- `adata.obsm['X_DR_sample']` — the same embedding as an ndarray.

This single key replaces the OLD two-key `X_DR_expression` / `X_DR_proportion` layout and the `(pseudo_dict, pseudo_adata)` return.

## Usage

```python
from sampledisco.sample_embedding import compute_sample_embedding

adata_integrated = compute_sample_embedding(
    adata_integrated,
    output_dir="/results/multiomics",
    sample_col="sample",
    celltype_col="cell_type",
    modality_col="modality",   # multi-omics: split units by modality
    pca_components=10,
    use_gpu=True,
)
```

!!! tip "GPU dispatch"
    Pass `use_gpu=True` (forwarded to the CPU/GPU dispatcher in `sampledisco.sample_embedding`) to use the RAPIDS-accelerated path when the RAPIDS stack is importable; otherwise SampleDisco falls back cleanly to CPU.
