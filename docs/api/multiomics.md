# Multi-omics API

The multi-omics branch builds one embedding of RNA and ATAC cells, derives the two views SampleDisco needs, types cells jointly, and lifts everything into one sample embedding (`uns['X_DR_sample']`) whose units are `<sample>_RNA` and `<sample>_ATAC`.

!!! tip "Run the wrapper, not the low-level functions"
    The supported entry point is the config-driven wrapper (`sampledisco -m complex --config <yaml>` with `run_multiomics_pipeline: true`). It runs the steps below in order and then [`downstream_analysis`](pipeline/downstream_analysis.md). The pages here are for reference and for step-by-step use.

| Function | Purpose |
| --- | --- |
| [`multiomics_preparation`](multiomics/multiomics_preparation.md) | scGLUE preprocessing and training, the embedding-only RNA+ATAC union (`preprocess/adata_sample.h5ad`), per-modality QC files for DGE, optional plots. Needs `scglue` and `bedtools`. |
| [`harmonize_xglue`](multiomics/harmonize_xglue.md) | Harmony on `obsm['X_glue']` with the sample (and batch) as covariates → the sample-removed `obsm['Z_comp']`; also aliases `X_glue` as `Z_rmd`. The wrapper runs it when the union has no composition key yet. |
| [`cell_types_multiomics`](multiomics/cell_types_multiomics_linux.md) | Leiden on the RNA cells in `Z_comp`, Jaccard-weighted SNN label transfer to the ATAC cells, then removal of clusters with almost no cells of one modality. |
| [`compute_sample_embedding`](multiomics/calculate_multiomics_sample_embedding.md) (multi-omics) | The shared sample embedding with `modality_col="modality"`. |
| [`run_autotune`](shared/run_autotune.md) | Picks the RMD weight α; `tune_on_modality` scores on one modality's units. Wrapper key: `multiomics_autotune_enable`. |
| [`save_embedding_to_h5ad`](shared/save_embedding_to_h5ad.md) | Stores `uns['X_DR_sample']` in the union file (the wrapper does this after the sample embedding). |
| [`visualize_multimodal_embedding`](downstream/visualize_multimodal_embedding.md) | Scatter plots of the `<sample>_<modality>` units. |
| [`replace_optimal_dimension_reduction`](multiomics/replace_optimal_dimension_reduction.md) | **Legacy.** Works only on the pre-0.3 two-key output layout; nothing in the current pipeline calls it. |

## Typical order (what the wrapper does)

```
multiomics_preparation          # scGLUE → X_glue; union preprocess/adata_sample.h5ad
   └─> harmonize_xglue          # Z_comp (skipped if a second scGLUE run already made it)
         └─> cell_types_multiomics
               └─> compute_sample_embedding(modality_col="modality")   # or run_autotune
                     └─> save_embedding_to_h5ad(union)
                           └─> downstream analyses (+ visualize_multimodal_embedding)
```

The embedding does not have to come from scGLUE. Any joint cell embedding works if you put the sample-preserved view in `obsm['Z_rmd']` and the sample-removed view in `obsm['Z_comp']`; the [MaxFuse tutorial](../tutorials/multiomics_maxfuse.md) builds both with [`harmonize_embedding`](utils/harmonize_embedding.md) and [`harmonize_xglue`](multiomics/harmonize_xglue.md).

For a runnable walkthrough, see the [Multi-omics tutorial](../tutorials/multiomics.md).
