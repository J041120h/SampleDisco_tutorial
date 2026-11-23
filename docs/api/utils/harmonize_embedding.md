# `harmonize_embedding`

Batch-corrects an `n_cells × n_dims` embedding with Harmony and returns the corrected array. It is the Harmony call used everywhere in SampleDisco's cell-level steps (both [preprocess](../rna/preprocess_linux.md) passes and [`harmonize_xglue`](../multiomics/harmonize_xglue.md)), and a convenient way to build `Z_rmd` / `Z_comp` for an embedding that SampleDisco did not produce (see the [MaxFuse tutorial](../../tutorials/multiomics_maxfuse.md)).

Backend: **harmony-pytorch** (`from harmony import harmonize`), a core dependency, with `n_jobs` = the CPUs usable by the process (its CPU-affinity mask, capped at `SLURM_CPUS_PER_TASK` when set). **harmonypy** is only a fallback for when harmony-pytorch cannot be imported; it runs with `device="cpu"` unless `use_gpu=True`. The sample-level correction inside [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) uses harmonypy directly and does not go through this function.

**Import:** `from sampledisco.utils.harmony_compat import harmonize_embedding`

**Source:** `utils/harmony_compat.py:25`

## Signature

```python
def harmonize_embedding(
    embedding: np.ndarray,
    obs,
    batch_key: Optional[Union[str, Sequence[str]]] = None,
    max_iter_harmony: int = 30,
    use_gpu: bool = False,
    seed: int = 42,
) -> np.ndarray
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `embedding` | array | — | `n_cells × n_dims` matrix, e.g. `adata.obsm["X_pca"]`. |
| `obs` | DataFrame | — | Cell metadata with the batch columns, row-aligned with `embedding` (usually `adata.obs`). |
| `batch_key` | str, list or None | `None` | Column(s) to correct for. `None` or an empty list returns the input unchanged (as float32). |
| `max_iter_harmony` | int | `30` | Maximum Harmony iterations. |
| `use_gpu` | bool | `False` | Run harmony-pytorch on the GPU. With the harmonypy fallback, `False` pins harmonypy ≥ 0.2 to the CPU; `True` lets it choose its device. |
| `seed` | int | `42` | Random seed. |

## Returns

`numpy.ndarray` of the same shape, dtype float32. The input is not modified; store the result yourself, e.g. `adata.obsm["Z_rmd"] = ...`.

## Usage

On the RNA chain's object, this repeats the second Harmony pass of [`preprocess`](../rna/preprocess_linux.md) (`X_pca_rmd` corrected for `batch`) twice:

```python
import anndata as ad
import numpy as np
from sampledisco.utils.harmony_compat import harmonize_embedding

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
Z = harmonize_embedding(adata.obsm["X_pca_rmd"], adata.obs, batch_key=["batch"], max_iter_harmony=30)
Z_again = harmonize_embedding(adata.obsm["X_pca_rmd"], adata.obs, batch_key=["batch"], max_iter_harmony=30)
print(Z.shape, Z.dtype)
print("max |Z - Z_again| =", float(np.abs(Z - Z_again).max()))
```

Output (harmony-pytorch 0.1.8; about 30 s per call on 8 CPUs, 1.9 GB):

```text
(29680, 20) float32
max |Z - Z_again| = 0.0
```

Repeated calls on the same machine give identical results. Running `preprocess` twice, in two separate processes on one node, also gave identical `Z_comp` and `Z_rmd`. Calling this function yourself on a reloaded object does not exactly reproduce the `Z_rmd` that `preprocess` computed inside its own run: in our runs the largest differences were 0.5 and 2.9.
