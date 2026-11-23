# ATAC API

The ATAC branch replaces RNA-style log-normalization with TF-IDF + LSI. Only preprocessing is ATAC-specific. Cell typing uses the same `cell_types` function as RNA: it detects `obsm['X_lsi']` and switches the neighbor graph to a cosine metric.

| Function | Import | Purpose |
| --- | --- | --- |
| [`preprocess`](atac/preprocess_linux.md) *(ATAC)* | `sampledisco.preparation.atac_preprocess_cpu` (GPU: `atac_preprocess_gpu.preprocess_gpu`) | QC, optional Scrublet doublet removal, TF-IDF + log1p, LSI (first component dropped by default), then Harmony. Writes `preprocess/adata_preprocessed.h5ad` with `obsm['Z_comp']` and `obsm['Z_rmd']`. |
| [`cell_types`](rna/cell_types_linux.md) *(shared)* | `sampledisco.preparation.cell_type_cpu` (GPU: `cell_type_gpu.cell_types_gpu`) | Leiden on `Z_comp`, cosine neighbors for ATAC. This is what the wrapper calls for ATAC. |

!!! warning "`cell_types_atac` is not the pipeline's ATAC typer"
    `sampledisco.preparation.ATAC_cell_type.cell_types_atac` is an older ATAC-only routine. The ATAC wrapper does not call it, and in 0.5.0 and 0.6.0 its `n_target_clusters` merge path calls a function that does not exist (`cell_type_dendrogram_atac`, a `NameError`). Use [`cell_types`](rna/cell_types_linux.md).

## Shared steps

- [`compute_sample_embedding`](shared/calculate_sample_embedding.md): the one sample-embedding entry point for RNA, ATAC and multi-omics (`from sampledisco.sample_embedding import compute_sample_embedding`). It writes `adata.uns['X_DR_sample']`; there is no `atac=` flag.
- [`run_autotune`](shared/run_autotune.md): picks the RMD weight α. In the wrapper, `atac_autotune_enable: true` replaces the plain sample-embedding call.
- [`build_sample_adata`](shared/build_sample_adata.md) and the [downstream](downstream.md) functions.

## Typical order

```
preprocess (ATAC)
   └─> cell_types                      # shared with RNA
          └─> compute_sample_embedding  (or run_autotune)
                 └─> build_sample_adata
                        └─> downstream analyses
```

The supported way to run this end to end is the config-driven wrapper with `run_atac_pipeline: true`:

```bash
sampledisco -m complex --config config.yaml
```

For a runnable walkthrough, see the [ATAC tutorial](../tutorials/atac.md).
