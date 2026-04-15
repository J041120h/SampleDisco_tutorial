# RNA API

The RNA branch of SampleDisc takes a cell-level `.h5ad` with raw counts and turns it into two sample-level embeddings. Two RNA-specific functions live in this section; sample embedding and resolution search are shared across modalities and live under **Shared**.

| Function | Purpose |
| --- | --- |
| [`preprocess_linux`](rna/preprocess_linux.md) | QC filter, normalize, log-transform, select HVGs, PCA, Harmony. Produces `adata_cell.h5ad` (Harmony-integrated for clustering) and `adata_sample.h5ad` (minimally processed for pseudobulk). |
| [`cell_types_linux`](rna/cell_types_linux.md) | Leiden clustering on the Harmony embedding, with optional recursive resolution tuning to hit a target cluster count and optional UMAP. Writes labels into `.obs["cell_type"]`. |

## Shared with ATAC

- [`calculate_sample_embedding`](shared/calculate_sample_embedding.md) — pseudobulk + dual DR.
- [`find_optimal_cell_resolution_linux`](shared/find_optimal_cell_resolution_linux.md) — CCA-guided resolution search.

## Typical order

```
preprocess_linux
   └─> cell_types_linux
          └─> calculate_sample_embedding
                 └─> downstream analyses (see Downstream section)
```

For a runnable walkthrough, see the [RNA tutorial](../tutorials/rna.md).
