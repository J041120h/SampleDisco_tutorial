# RNA API

The RNA branch of SampleDisco takes a cell-level `.h5ad` with raw counts and turns it into a single sample-level embedding (`uns['X_DR_sample']`). Two RNA-specific functions live in this section; sample embedding is shared across modalities and lives under **Shared**.

| Function | Purpose |
| --- | --- |
| [`preprocess`](rna/preprocess_linux.md) | QC filter, normalize, log-transform, select HVGs, PCA, 2-pass Harmony. Writes one `adata_preprocessed.h5ad` carrying `obsm['Z_clust']` (sample-removed) and `obsm['Z_rmd']` (sample-preserved). |
| [`cell_types`](rna/cell_types_linux.md) | Leiden clustering on the sample-removed `Z_clust` embedding, with optional adaptive resolution tuning to hit a target cluster count and optional UMAP. Writes labels into `.obs["cell_type"]`. |

## Shared with ATAC

- [`compute_sample_embedding`](shared/calculate_sample_embedding.md) — composition blocks on `Z_clust` + RMD displacement block on `Z_rmd`, reduced to the single key `uns['X_DR_sample']`.
- [`run_autotune`](shared/run_autotune.md) — parameter selection: optimizes the RMD-weight α (composition vs displacement), enabled via the `rna_autotune_enable` flag in the config-driven wrapper.

## Typical order

```
preprocess
   └─> cell_types
          └─> compute_sample_embedding
                 └─> downstream analyses (see Downstream section)
```

For a runnable walkthrough, see the [RNA tutorial](../tutorials/rna.md).
