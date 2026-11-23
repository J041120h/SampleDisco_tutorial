# `harmonize_embedding`

Batch-corrects an `n_cells × n_dims` embedding with Harmony and returns the corrected array. It is the Harmony call used everywhere in SampleDisco's cell-level steps (both [preprocess](../rna/preprocess_linux.md) passes and [`harmonize_xglue`](../multiomics/harmonize_xglue.md)), and a convenient way to build `Z_rmd` / `Z_comp` for an embedding that SampleDisco did not produce (see the [MaxFuse tutorial](../../tutorials/multiomics_maxfuse.md)).

Backend: **harmony-pytorch** (`from harmony import harmonize`) when it and PyTorch import, otherwise **harmonypy**. In 0.6.0, harmony-pytorch gets `n_jobs` = the CPUs usable by the process (its CPU-affinity mask, capped at `SLURM_CPUS_PER_TASK` when set) instead of every core of the node, which also makes it deterministic. harmonypy ≥ 0.2 has a PyTorch backend that would otherwise pick a visible GPU, so it gets `device="cpu"` unless `use_gpu=True`. The sample-level correction inside [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) uses harmonypy directly and does not go through this function.

**Import:** `from sampledisco.utils.harmony_compat import harmonize_embedding`

**Source:** `utils/harmony_compat.py:26`

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

On the RNA chain's object, this repeats the second Harmony pass of [`preprocess`](../rna/preprocess_linux.md) (`X_pca_rmd` corrected for `batch`) and compares it with the stored `Z_rmd`:

```python
import anndata as ad
import numpy as np
from sampledisco.utils.harmony_compat import harmonize_embedding

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
Z = harmonize_embedding(adata.obsm["X_pca_rmd"], adata.obs, batch_key=["batch"], max_iter_harmony=30)
print(Z.shape, Z.dtype)
print("max |Z - Z_rmd| =", float(np.abs(Z - adata.obsm["Z_rmd"]).max()))
```

Output (about 7 s; this environment has no harmony-pytorch, so harmonypy 0.2.0 ran, on the CPU):

```text
(29680, 20) float32
max |Z - Z_rmd| = 0.14886808395385742
```

The result is close to, but not bit-identical with, the `Z_rmd` that `preprocess` stored in an earlier job. Across three runs the largest difference was 0.011, 0.149 and 0.149, so harmonypy results are not bit-reproducible between runs. The deterministic setting in 0.6.0 concerns harmony-pytorch, which this environment does not have.
