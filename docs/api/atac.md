# ATAC API

The ATAC branch replaces RNA-style log-normalization with TF-IDF + LSI. Only one function is ATAC-specific; cell typing and sample embedding are shared with RNA.

| Function | Purpose |
| --- | --- |
| [`preprocess_linux`](atac/preprocess_linux.md) *(ATAC)* | TF-IDF normalization, LSI projection, optional doublet detection, Harmony on LSI. Drops the first LSI component (depth-correlated) by default. |

## Shared with RNA

- [`cell_types_linux`](rna/cell_types_linux.md) — auto-detects `X_lsi_harmony` in `.obsm` and switches to cosine-metric neighborhoods when it sees an ATAC embedding.
- [`calculate_sample_embedding`](shared/calculate_sample_embedding.md) — set `atac=True` to enable TF-IDF/LSI-aware pseudobulk aggregation.
- [`find_optimal_cell_resolution_linux`](shared/find_optimal_cell_resolution_linux.md) — set `modality="atac"`.

## Typical order

```
preprocess_linux (ATAC)
   └─> cell_types_linux               # same function as RNA
          └─> calculate_sample_embedding(atac=True)
                 └─> downstream analyses
```

For a runnable walkthrough, see the [ATAC tutorial](../tutorials/atac.md).
