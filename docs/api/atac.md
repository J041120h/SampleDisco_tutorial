# ATAC API

The ATAC branch replaces RNA-style log-normalization with TF-IDF + LSI. Preprocessing and cell typing are ATAC-specific; the sample embedding is shared with RNA.

| Function | Purpose |
| --- | --- |
| [`preprocess`](atac/preprocess_linux.md) *(ATAC)* | TF-IDF normalization, LSI projection, optional doublet detection, Harmony on LSI. Drops the first LSI component (depth-correlated) by default. Imported as `from sampledisco.preparation.atac_preprocess_cpu import preprocess`. |
| [`cell_types_atac`](rna/cell_types_linux.md) *(ATAC)* | Leiden cell typing on the ATAC DR embedding, with dendrogram / differential-peak analysis. Imported as `from sampledisco.preparation.ATAC_cell_type import cell_types_atac`. |

## Shared with RNA

- [`compute_sample_embedding`](shared/calculate_sample_embedding.md) — the single, unified sample-embedding entry point for RNA, ATAC, and multi-omics (`from sampledisco.sample_embedding import compute_sample_embedding`). It writes the final sample embedding to `adata.uns['X_DR_sample']`; there is no longer an `atac=True` flag.
- [`run_autotune`](shared/run_autotune.md) — parameter selection for the sample embedding (`from sampledisco.parameter_selection.autotune import run_autotune`). It optimizes the RMD-weight alpha (composition vs displacement), enabled via the `atac_autotune_enable` flag in the config-driven wrapper.

## Typical order

```
preprocess (ATAC)
   └─> cell_types_atac                # ATAC-specific cell typing
          └─> compute_sample_embedding
                 └─> downstream analyses
```

The supported way to run this end-to-end is the config-driven wrapper:

```bash
sampledisco -m complex --config config.yaml
```

For a runnable walkthrough, see the [ATAC tutorial](../tutorials/atac.md).
