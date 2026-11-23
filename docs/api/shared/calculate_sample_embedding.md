# `compute_sample_embedding`

Lifts a cell-level embedding into one sample-level embedding. For each **unit** (a sample; in multi-omics a `<sample>_<modality>` pair) it builds four blocks:

1. **Composition blocks**, on the sample-removed embedding `Z_comp`:
    - **A1**: cell-type composition (one-hot `celltype_col`, mean per unit).
    - **A2**: soft k-means composition with `K_med = min(medium_K, max(2, n_cells // 200))` centers.
    - **A3**: soft k-means composition with `K_fine = min(fine_K, max(2, n_cells // 100))` centers.
2. **RMD displacement block** (`use_rmd=True`), on the sample-preserved embedding `Z_rmd`: for each cell type, the mean of the unit's cells minus the mean of that cell type's cells in the *other* units of the same group (leave-one-out), down-weighted for units with few cells of that type, then reduced by PCA to at most `rmd_dim_per_cluster` dimensions per cell type. The group is the modality in multi-omics, the unit's majority `batch_col` label in single-omics, or one group for all units when there is no `batch_col`.

Each block is centered and scaled to Frobenius norm √N · weight, the blocks are concatenated, and PCA keeps `min(pca_components, n_units − 1)` components. A sample-level batch correction follows (see `batch_method`). The result is stored in `adata.uns['X_DR_sample']`, the single key every downstream function reads.

The same function serves RNA, ATAC and multi-omics; pass `modality_col="modality"` for multi-omics (see the [multi-omics page](../multiomics/calculate_multiomics_sample_embedding.md)). The tutorial [Sample embedding in depth](../../tutorials/sample_embedding.md) walks through the blocks, weights and batch options on the demo.

**Import:** `from sampledisco.sample_embedding import compute_sample_embedding` (also `sampledisco.compute_sample_embedding`).

**Source:** `sample_embedding/__init__.py:23` (CPU/GPU dispatcher) → `sample_embedding/sample_embedding.py:143` (CPU) or `sample_embedding/sample_embedding_gpu.py:162` (GPU).

!!! note "Renamed and unified"
    This replaces the old `calculate_sample_embedding` / `calculate_multiomics_sample_embedding` pair, which built two keys (`X_DR_expression`, `X_DR_proportion`) and returned `(pseudo_dict, pseudo_adata)`. The current function writes one key and returns the cell-level AnnData.

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

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level object with `obs[sample_col]`, `obs[celltype_col]` and the two embeddings. Modified in place. |
| `output_dir` | str | — | Root output directory (see Output files). |
| `use_gpu` | bool | `False` | Use the RAPIDS implementation: the soft assignments and the final PCA run on the GPU (`cupy`/`cuml`). The k-means centres are fitted on the CPU with scikit-learn's `MiniBatchKMeans` (cuML's if the build has one), and the sample-level Harmony runs on the CPU with harmonypy, as on the CPU path. Falls back to CPU with a message when `cupy`/`cuml` are missing or the GPU call fails. |
| `sample_col` | str | `"sample"` | Sample column. |
| `celltype_col` | str | `"cell_type"` | Cell-type column for A1 and the RMD block. At least 2 cell types and 2 units are required. |
| `comp_emb_key` | str or None | `None` | `.obsm` key of the sample-removed view. `None` resolves `Z_comp`, then `Z_clust`, `X_pca_harmony`, `X_lsi_harmony`, `X_glue_harmony` (legacy names give a `FutureWarning`). An explicit `"Z_comp"` on an old file falls back the same way in 0.6.0 (with a `FutureWarning`); 0.5.0 raised `KeyError`. See [key resolvers](../utils/embedding_keys.md). |
| `rmd_emb_key` | str or None | `None` | `.obsm` key of the sample-preserved view. `None` resolves `Z_rmd`, then `Z_cmd`, `X_pca_harmony_nosamp`, `X_lsi_harmony_nosamp`. If none exists, it falls back to the composition key with a `RuntimeWarning`; the RMD block is then near-degenerate and the result is not the published method. There is no `X_glue` fallback here (the multi-omics wrapper adds one). |
| `modality_col` | str or None | `None` | Set (e.g. `"modality"`) for multi-omics: units become `<sample>_<modality>` and the RMD reference group is the modality. |
| `batch_col` | str, list or None | `None` | Sample-level batch column(s). The first one sets the single-omics RMD group (majority label per unit) and the Harmony labels; with two or more columns, Harmony gets one covariate per column (majority label per unit). |
| `medium_K` | int | `120` | Upper bound for `K_med` (A2). |
| `fine_K` | int | `300` | Upper bound for `K_fine` (A3). |
| `rmd_dim_per_cluster` | int | `8` | Maximum RMD dimensions per cell type. |
| `use_clr` | bool | `False` | CLR-transform the three composition blocks. |
| `use_rmd` | bool | `True` | Include the RMD block. `False` gives a composition-only embedding (3 blocks). |
| `block_weights` | list or None | `None` | Explicit weights, one per block (4 with RMD, 3 without). `None` uses w_A1 = √(K_fine / K_c), w_A2 = √(K_fine / K_med), w_A3 = 1, w_RMD = α, where K_c is the number of cell types. |
| `rmd_weight` | float or `"equal"` | `"equal"` | α. `"equal"` gives the RMD block the same energy as the three composition blocks together: α² = w_A1² + w_A2² + w_A3², rounded to 2 decimals. A finite positive number fixes α (0.60 was the default before 0.5.0). Anything else raises `ValueError`. [`run_autotune`](run_autotune.md) searches α instead. |
| `pca_components` | int | `10` | PCs kept; capped at `n_units − 1`. |
| `batch_method` | str | `"harmony"` | Sample-level correction: `"harmony"` (harmonypy; falls back to linear regression, then to no correction, if it fails), `"linear"` (regress out the batch labels per PC) or `"none"`. Any correction runs only when there are **at least 8 units and at least 2 batch labels**; otherwise the PCA scores are used as is. Without `batch_col`, single-omics units all share one label, so no correction runs; multi-omics units are labelled by modality (or modality × batch). |
| `save` | bool | `True` | Write the CSV and block files (and update the h5ad, see `save_cell_adata`). |
| `verbose` | bool | `True` | Print progress. |
| `seed` | int | `42` | Seed for k-means, RMD, PCA and Harmony. |
| `cluster_emb_key` | str or None | `None` | Deprecated 0.2.0 spelling of `comp_emb_key` (`FutureWarning`); used only when `comp_emb_key` is `None`. |
| `save_cell_adata` | bool | `True` | With `save=True`, also store the two `uns` entries in `{output_dir}/preprocess/adata_preprocessed.h5ad` **if that file exists**, through [`save_embedding_to_h5ad`](save_embedding_to_h5ad.md): only those entries are replaced when the file holds the same object as memory, otherwise the whole file is rewritten. `False` skips this and warns if the file still holds an older `uns['X_DR_sample']`. Config key: `*_save_cell_adata_after_embedding`. |

## Returns

The **same cell-level `AnnData`**, modified in place, with:

- `.uns['X_DR_sample']`: `pandas.DataFrame`, units × PCs (`PC1`, `PC2`, …), index named `sample`.
- `.uns['sample_embedding_params']`: dict with `medium_K`, `fine_K` (the K values actually used), `K_c`, `use_clr`, `use_rmd`, `rmd_weight` (the resolved α), `block_weights`, `pca_components`, `batch_method`, `comp_emb_key`, `rmd_emb_key`, `modality_col`, `batch_col`, `batch_cols_multi`, `seed` and `backend` (`"cpu"` or `"gpu"`).

There is no `.obsm['X_DR_sample']` on the cell-level object. [`build_sample_adata`](build_sample_adata.md) turns the result into one row per unit for the downstream functions.

## Output files

With `save=True`:

- `{output_dir}/sample_embedding/sample_embedding.csv`: the embedding.
- `{output_dir}/sample_embedding/sample_embedding_blocks.npz`: `unit_ids`, `A1_cell_types`, `A1`, `A2`, `A3`, `RMD` (the unweighted blocks).
- `{output_dir}/preprocess/adata_preprocessed.h5ad`: updated only if it already exists and `save_cell_adata=True`.

## Usage

Continuing from the [`cell_types`](../rna/cell_types_linux.md) example. The demo's two processing batches (`Aruna`, `Wilk`) and 8 samples meet the Harmony condition:

```python
import anndata as ad
from sampledisco.sample_embedding import compute_sample_embedding

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
adata = compute_sample_embedding(
    adata,
    "api_demo/rna",
    sample_col="sample",
    celltype_col="cell_type",
    batch_col="batch",
)
params = adata.uns["sample_embedding_params"]
print(adata.uns["X_DR_sample"].shape)
print(params["K_c"], params["medium_K"], params["fine_K"], params["rmd_weight"])
print([round(w, 2) for w in params["block_weights"]])
```

**Writes** → `api_demo/rna/sample_embedding/sample_embedding.csv`, `sample_embedding_blocks.npz`, and `uns['X_DR_sample']` into `api_demo/rna/preprocess/adata_preprocessed.h5ad` (updated in place). Printed at the end:

```text
(8, 7)
8 120 296 6.36
[6.08, 1.57, 1.0, 6.36]
```

8 cell types; `K_fine` is capped at 29,680 cells // 100 = 296; `"equal"` resolves to α = 6.36 = √(6.08² + 1.57² + 1²); the embedding has 7 PCs, not 10, because 8 units allow at most 7. The log line `[batch correction] 2 groups (group only)` shows that sample-level Harmony ran on the two batches. About 3 s and 2.1 GB.

For multi-omics, pass `modality_col="modality"`; see [`compute_sample_embedding` (multi-omics)](../multiomics/calculate_multiomics_sample_embedding.md).

!!! tip "Prefer the config-driven wrapper"
    The wrapper calls this function with `*_sample_embedding_medium_K`, `*_sample_embedding_fine_K`, `*_sample_embedding_rmd_dim`, `*_sample_embedding_use_clr`, `*_sample_embedding_use_rmd`, `*_sample_embedding_block_weights`, `*_sample_embedding_rmd_weight`, `*_sample_embedding_pca_components`, `*_sample_embedding_batch_method`, with `batch_col` from `rna_sample_level_batch_col` / `atac_sample_level_batch_col` / `multiomics_batch_col`. GPU is used only when `use_gpu: true`.
