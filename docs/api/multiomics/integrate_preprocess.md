# `integrate_preprocess`

!!! warning "Removed"
    `integrate_preprocess` no longer exists as a standalone function in the current `sampledisco` package. The OLD per-modality QC + union-build step that read and overwrote `{output_dir}/preprocess/adata_sample.h5ad` has been **folded into [`multiomics_preparation`](multiomics_preparation.md)**. There is no drop-in replacement to call directly — use the config-driven wrapper or `multiomics_preparation` instead.

## What replaced it

The post-integration QC and embedding-union assembly that this function used to perform are now handled inside `multiomics_preparation`, controlled by two flags:

| Flag | Default | Role |
| --- | --- | --- |
| `run_preprocess_per_modality` | `True` | Per-modality QC + normalize for downstream (the work `integrate_preprocess` used to do). |
| `run_merge` | `True` | Build the embedding-only union AnnData across RNA and ATAC. |

Internally these steps call the helpers in `sampledisco.preparation.multi_omics_merge`:

- `build_embedding_union` — assembles the merged, embedding-carrying AnnData.
- `preprocess_rna_for_downstream` — per-modality RNA QC + normalize.
- `preprocess_atac_for_downstream` — per-modality ATAC QC + normalize.

## Recommended approach

Run the config-driven wrapper, which orchestrates the full multi-omics pipeline (preprocess → scGLUE integration → cell typing → sample embedding → downstream):

```bash
sampledisco -m complex --config config.yaml
```

Or call `multiomics_preparation` directly, leaving `run_preprocess_per_modality` and `run_merge` enabled:

```python
from sampledisco.preparation.multi_omics_glue import multiomics_preparation

multiomics_preparation(
    rna_file="/data/rna.h5ad",
    atac_file="/data/atac.h5ad",
    run_preprocess_per_modality=True,
    run_merge=True,
)
```

See [`multiomics_preparation`](multiomics_preparation.md) for the full parameter list.
