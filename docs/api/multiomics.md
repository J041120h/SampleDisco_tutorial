# Multi-omics API

The multi-omics branch integrates scRNA + scATAC via scGLUE, cleans the merged object, assigns joint cell types, and lifts the joint cell embedding into a single sample embedding (`uns['X_DR_sample']`) shared with the single-modality pipelines.

!!! tip "Run the wrapper, not the low-level functions"
    The supported entry point is the config-driven wrapper (`sampledisco -m complex --config <yaml>` with `run_multiomics_pipeline=True`), which orchestrates preparation → joint cell typing → `compute_sample_embedding` → downstream analysis for you. The functions below are documented for reference and for advanced step-by-step use.

| Function | Purpose |
| --- | --- |
| [`multiomics_preparation`](multiomics/multiomics_preparation.md) | Full scGLUE pipeline — preprocessing, adversarial training, per-modality QC, embedding-union merge, optional visualization. |
| [`cell_types_multiomics`](multiomics/cell_types_multiomics_linux.md) | Leiden on RNA cells in the joint embedding, then Jaccard-weighted SNN label transfer to ATAC. |
| [`compute_sample_embedding`](multiomics/calculate_multiomics_sample_embedding.md) | Single-key sample embedding for multi-omics (`uns['X_DR_sample']`); pass `modality_col='modality'`. |
| [`replace_optimal_dimension_reduction`](multiomics/replace_optimal_dimension_reduction.md) | Legacy optimal-resolution merge helper. |
| [`run_autotune`](shared/run_autotune.md) | Parameter selection — optimizes the RMD-weight α (composition vs displacement); pass `modality_col='modality'`. Enabled via `multiomics_autotune_enable` in the wrapper. |

## Typical order

```
multiomics_preparation                  # scGLUE integration + per-modality QC + merge
   └─> cell_types_multiomics            # joint cell typing
         └─> compute_sample_embedding   # writes uns['X_DR_sample'], modality_col='modality'
               └─> (optional) multiomics_autotune_enable   # alpha / block-weight tuning
                     └─> downstream analyses
```

For a runnable walkthrough, see the [Multi-omics tutorial](../tutorials/multiomics.md).
