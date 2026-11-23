# RNA API

The RNA branch takes a cell-level `.h5ad` with raw counts and turns it into one sample embedding, `uns['X_DR_sample']`. Two functions are RNA-specific in their inputs; everything from the sample embedding on is shared with ATAC and multi-omics.

| Function | Import | Purpose |
| --- | --- | --- |
| [`preprocess`](rna/preprocess_linux.md) | `sampledisco.preparation.rna_preprocess_cpu` (GPU: `rna_preprocess_gpu.preprocess_gpu`) | QC, normalize + log1p, HVGs, PCA, then Harmony. Writes `preprocess/adata_preprocessed.h5ad` with `obsm['Z_comp']` (sample removed) and `obsm['Z_rmd']` (sample preserved). |
| [`cell_types`](rna/cell_types_linux.md) | `sampledisco.preparation.cell_type_cpu` (GPU: `cell_type_gpu.cell_types_gpu`) | Leiden on `Z_comp`, optionally merged down to a target number of clusters. Labels go to `obs['cell_type']`. The same function types ATAC cells. |

## Shared steps

- [`compute_sample_embedding`](shared/calculate_sample_embedding.md): composition blocks on `Z_comp` + RMD displacement block on `Z_rmd` → `uns['X_DR_sample']`.
- [`run_autotune`](shared/run_autotune.md): picks the RMD weight α instead of using the default. In the wrapper, `rna_autotune_enable: true` replaces the plain sample-embedding call.
- [`build_sample_adata`](shared/build_sample_adata.md): one-row-per-sample AnnData for the downstream functions.
- [Downstream](downstream.md) functions and [`downstream_analysis`](pipeline/downstream_analysis.md).

## Typical order

```
preprocess
   └─> cell_types
          └─> compute_sample_embedding   (or run_autotune)
                 └─> build_sample_adata
                        └─> downstream analyses
```

The wrapper runs the same order from a config (`run_rna_pipeline: true`); GPU versions are used only when the config sets `use_gpu: true` and RAPIDS imports. For a runnable walkthrough, see the [RNA tutorial](../tutorials/rna.md).
