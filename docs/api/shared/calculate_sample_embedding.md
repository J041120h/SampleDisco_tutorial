# `compute_sample_embedding`

!!! note "Renamed and unified"
    This is the current core sample-embedding function. It replaces the old
    `calculate_sample_embedding` / `calculate_multiomics_sample_embedding`
    pair. The old API built a **two-key** embedding (`X_DR_expression` +
    `X_DR_proportion`) and returned `(pseudo_dict, pseudo_adata)`. The current
    API builds a **single-key** embedding written in place to
    `adata.uns['X_DR_sample']` and **returns the modified cell-level AnnData**.

Lifts a cell-level embedding into a single sample-level embedding (the SampleDisco "singleRMD" recipe). For each unit (a sample, or a `<sample>_RNA` / `<sample>_ATAC` unit in multi-omics) it combines:

1. **Multi-resolution composition blocks** computed on the sample-removed embedding `Z_clust`:
   - **A1** — coarse cell-type composition (one-hot, mean per unit).
   - **A2** — soft k-means composition at `medium_K`.
   - **A3** — soft k-means composition at `fine_K`.
2. **RMD displacement block** (`use_rmd=True`) computed on the sample-preserved embedding `Z_rmd` — each coarse cell type's mean position relative to a leave-one-out reference built from comparable units. Captures within-cell-type state shifts.

The blocks are inverse-variance weighted, Frobenius-stacked, PCA-reduced (`pca_components`), and Harmony-corrected at the sample level. The result is written to `adata.uns['X_DR_sample']` (a pandas DataFrame, units × PCs) and consumed by every downstream module.

Works for RNA, ATAC, and multi-omics: pass `modality_col='modality'` for the multi-omics case — there is no separate multi-omics entry point. GPU acceleration is enabled via `use_gpu=True`, which dispatches to the RAPIDS implementation.

**Source:** `sample_embedding/__init__.py` — `compute_sample_embedding`, the `use_gpu` dispatcher. The CPU and GPU implementations it delegates to live in `sample_embedding/sample_embedding.py` and `sample_embedding/sample_embedding_gpu.py` respectively.

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
| `adata` | AnnData | — | Cell-level AnnData with the cluster embedding in `.obsm[cluster_emb_key]` and a cell-type column in `.obs`. Mutated in place. |
| `output_dir` | str | — | Writes the embedding CSV to `{output_dir}/sample_embedding/` and re-saves `{output_dir}/preprocess/adata_preprocessed.h5ad`. |
| `use_gpu` | bool | `False` | Dispatch to the RAPIDS GPU implementation. |
| `sample_col` | str | `"sample"` | Column in `.obs` identifying samples. |
| `celltype_col` | str | `"cell_type"` | Column in `.obs` identifying cell types. |
| `cluster_emb_key` | str | `"Z_clust"` | obsm key for the **sample-removed** embedding (composition blocks). |
| `rmd_emb_key` | str, optional | `None` | obsm key for the **sample-preserved** embedding (RMD block). When `None`, resolves to `"Z_rmd"` if present, else falls back to `cluster_emb_key`. |
| `modality_col` | str, optional | `None` | Set to `"modality"` for multi-omics; defines the RMD grouping and the `<sample>_RNA` / `<sample>_ATAC` unit ids. |
| `batch_col` | str or list, optional | `None` | Sample-level batch column(s). The first labels the RMD grouping; supplying ≥2 enables multi-covariate sample-level Harmony. |
| `medium_K` | int | `120` | Target k for the A2 medium-resolution soft k-means (capped by cell count). |
| `fine_K` | int | `300` | Target k for the A3 fine-resolution soft k-means (capped by cell count). |
| `rmd_dim_per_cluster` | int | `8` | Max RMD dimensions retained per coarse cluster. |
| `use_clr` | bool | `False` | CLR-transform the composition blocks. |
| `use_rmd` | bool | `True` | Include the RMD displacement block. |
| `block_weights` | list, optional | `None` | Explicit per-block weights. When `None`, weights are auto-derived from `K_c`/`medium_K`/`fine_K` via an inverse-variance schedule. |
| `rmd_weight` | float | `0.60` | Relative weight of the RMD block in the auto-derived schedule. |
| `pca_components` | int | `10` | Number of PCA components in the final sample embedding. |
| `batch_method` | str | `"harmony"` | Sample-level batch-correction method applied after PCA. |
| `save` | bool | `True` | Write the embedding CSV and re-save the preprocessed h5ad. |
| `verbose` | bool | `True` | Print progress. |
| `seed` | int | `42` | Random seed for k-means / RMD / Harmony. |

## Returns

`AnnData` — the **same cell-level AnnData**, mutated in place. It carries:

- `.uns['X_DR_sample']` — a pandas DataFrame (units × `pca_components`); the single sample embedding consumed by all downstream modules.
- `.uns['sample_embedding_params']` — a dict of the parameters and derived cluster counts (`K_c`, `medium_K`, `fine_K`, `block_weights`, etc.).

!!! tip "Materializing a per-sample AnnData"
    Downstream consumers that expect one row per sample can build it from the
    embedding with `build_sample_adata`:

    ```python
    from sampledisco.sample_embedding.sample_embedding import build_sample_adata

    pseudo_adata = build_sample_adata(adata, sample_col="sample")
    ```

## Output files

- `{output_dir}/sample_embedding/sample_embedding.csv`
- `{output_dir}/preprocess/adata_preprocessed.h5ad` (re-saved with `.uns['X_DR_sample']`, if it already exists)

## Usage

```python
from sampledisco.sample_embedding import compute_sample_embedding

adata = compute_sample_embedding(
    adata,
    output_dir="sampledisco_demo_output/rna",
    sample_col="sample",
    celltype_col="cell_type",
    cluster_emb_key="Z_clust",
    rmd_emb_key="Z_rmd",
    pca_components=10,
    use_gpu=True,
)

sample_embedding = adata.uns["X_DR_sample"]  # DataFrame, units × PCs
```

For multi-omics, pass `modality_col='modality'`:

```python
adata = compute_sample_embedding(
    adata,
    output_dir="sampledisco_demo_output/multiomics",
    modality_col="modality",
    cluster_emb_key="Z_clust",
    rmd_emb_key="Z_rmd",
)
```
