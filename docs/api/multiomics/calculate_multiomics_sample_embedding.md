# `compute_sample_embedding` (multi-omics)

There is no separate multi-omics function: the shared [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) handles RNA, ATAC and multi-omics. This page covers what changes when you pass `modality_col="modality"`. The full parameter table, the block definitions and the output files are on the [shared page](../shared/calculate_sample_embedding.md).

**Import:** `from sampledisco.sample_embedding import compute_sample_embedding`

**Source:** `sample_embedding/__init__.py:23` (CPU/GPU dispatcher) → `sample_embedding/sample_embedding.py:162` (CPU) or `sample_embedding/sample_embedding_gpu.py:164` (GPU).

!!! warning "Renamed and unified"
    The old `calculate_multiomics_sample_embedding` and its two-key output (`X_DR_expression`, `X_DR_proportion`) are gone. The current function writes one key, `uns['X_DR_sample']`.

## Signature

```python
def compute_sample_embedding(
    adata: AnnData,
    output_dir: str,
    *,
    use_gpu: bool = False,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    comp_emb_key: Optional[str] = None,
    rmd_emb_key: Optional[str] = None,
    modality_col: Optional[str] = None,
    batch_col: Optional[Union[str, List[str]]] = None,
    medium_K: int = 120,
    fine_K: int = 300,
    rmd_dim_per_cluster: int = 8,
    use_clr: bool = False,
    use_rmd: bool = True,
    block_weights: Optional[List[float]] = None,
    rmd_weight: Union[float, str] = "equal",
    pca_components: int = 10,
    batch_method: str = "harmony",
    save: bool = True,
    verbose: bool = True,
    seed: int = 42,
    cluster_emb_key: Optional[str] = None,
    save_cell_adata: bool = True,
) -> AnnData
```

## What changes with `modality_col`

| Aspect | Multi-omics behaviour |
| --- | --- |
| Units | One per (sample, modality) pair present, with id `<sample>_<modality>`. If the sample ids already end in `_<modality>` (as in the demo file), the suffix is not doubled. |
| RMD reference group | The modality: each RNA unit is compared with the other RNA units, each ATAC unit with the other ATAC units. |
| Sample-level correction labels | The modality, or `<modality>__<batch>` when `batch_col` is given (and that splits the units into fewer groups than there are units). With two or more `batch_col` columns, one covariate per column instead. Correction needs at least 8 units. |
| Embedding keys | `comp_emb_key=None` resolves `Z_comp` (or a legacy name). `rmd_emb_key=None` resolves `Z_rmd`, then `Z_cmd` and the other legacy names; **there is no `X_glue` fallback inside this function**. If no sample-preserved key is found, it falls back to the composition key with a `RuntimeWarning`. The multi-omics wrapper resolves the keys itself and does fall back to `X_glue`; [`harmonize_xglue`](harmonize_xglue.md) also aliases `X_glue` as `Z_rmd`. |
| `celltype_col` | Must cover both modalities, e.g. from [`cell_types_multiomics`](cell_types_multiomics_linux.md). |

## Returns

The same cell-level AnnData, modified in place, with `uns['X_DR_sample']` (a `DataFrame`, units × PCs, units `<sample>_RNA` / `<sample>_ATAC`) and `uns['sample_embedding_params']`. It does **not** set `obsm['X_DR_sample']`; [`build_sample_adata`](../shared/build_sample_adata.md) with `modality_col="modality"` gives the per-unit AnnData that has it.

## Output files

The same as on the [shared page](../shared/calculate_sample_embedding.md#output-files): `{output_dir}/sample_embedding/sample_embedding.csv` and `sample_embedding_blocks.npz`. The automatic h5ad update targets `{output_dir}/preprocess/adata_preprocessed.h5ad`, which does not exist in a multi-omics run; store the embedding in the union file with [`save_embedding_to_h5ad`](../shared/save_embedding_to_h5ad.md).

## Usage

Multi-omics chain, after [`cell_types_multiomics`](cell_types_multiomics_linux.md):

```python
import anndata as ad
from sampledisco.sample_embedding import compute_sample_embedding
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

mo = ad.read_h5ad("api_demo/multiomics/preprocess/adata_sample.h5ad")
mo = compute_sample_embedding(
    mo,
    "api_demo/multiomics",
    sample_col="sample",
    celltype_col="cell_type",
    modality_col="modality",   # units become <sample>_RNA / <sample>_ATAC
)
pseudo_mo = build_sample_adata(mo, sample_col="sample", modality_col="modality")
print(mo.uns["X_DR_sample"].shape)
print(pseudo_mo.obs[["modality", "sev.level"]].head())
```

**Writes** → `api_demo/multiomics/sample_embedding/sample_embedding.csv` and `sample_embedding_blocks.npz`. From the log and the prints:

```text
[sample_embedding] 16 units; 59010 cells; comp_emb dim=50
[A1] coarse cell-type composition: shape=(16, 16)
...
[sample_embedding] weights=[4.33, 1.581, 1.0, 4.72] (n_blocks=4)
  [batch correction] 2 groups (group only)
...
(16, 10)
                 modality sev.level
sample
CoV-3-Wilk_RNA        RNA         4
CoV-34-Aruna_RNA      RNA         4
CoV-36-Aruna_RNA      RNA         4
CoV-4-Wilk_RNA        RNA         4
HD-1-Wilk_RNA         RNA         1
```

16 units (8 RNA + 8 ATAC samples; the demo's sample ids already end in `_RNA` / `_ATAC`, so they are not doubled), 16 cell types, α = 4.72 from `"equal"`, 10 PCs (16 units allow up to 15), and sample-level Harmony on the two modalities. About 6 s and 1.8 GB. Nothing is written to the union file; see [`save_embedding_to_h5ad`](../shared/save_embedding_to_h5ad.md).

!!! tip "GPU"
    `use_gpu=True` uses the RAPIDS implementation when `cupy` and `cuml` import, and otherwise prints a message and runs on CPU. On the GPU path, the soft assignments and the final PCA run on the GPU; the k-means centres and the sample-level Harmony run on the CPU (see the [shared page](../shared/calculate_sample_embedding.md)). In the wrapper, the multi-omics branch uses the GPU only when both `use_gpu: true` and `multiomics_use_gpu: true` (the latter is the default).
