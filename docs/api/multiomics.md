# Multi-omics API

The multi-omics branch integrates scRNA + scATAC via GLUE, cleans the merged object, assigns joint cell types, and builds sample embeddings in the same dual-embedding format as the single-modality pipelines.

| Function | Purpose |
| --- | --- |
| [`multiomics_preparation`](multiomics/multiomics_preparation.md) | Full GLUE pipeline — preprocessing, adversarial training, gene-activity scoring, optional visualization. |
| [`integrate_preprocess`](multiomics/integrate_preprocess.md) | Post-GLUE QC filtering and per-modality metadata merge. |
| [`cell_types_multiomics_linux`](multiomics/cell_types_multiomics_linux.md) | Leiden on RNA cells in the joint embedding, then cosine k-NN label transfer to ATAC. |
| [`calculate_multiomics_sample_embedding`](multiomics/calculate_multiomics_sample_embedding.md) | Multi-omics pseudobulk + dual DR (`X_DR_expression`, `X_DR_proportion`). |
| [`find_optimal_cell_resolution_multiomics_linux`](multiomics/find_optimal_cell_resolution_multiomics_linux.md) | CCA-guided search for the Leiden resolution that best aligns with a phenotype. |
| [`replace_optimal_dimension_reduction`](multiomics/replace_optimal_dimension_reduction.md) | Merges the optimal embeddings back into `pseudobulk_sample.h5ad`. |

## Typical order

```
multiomics_preparation                      # GLUE integration
   └─> integrate_preprocess                 # QC on the merged object
         └─> cell_types_multiomics_linux    # joint cell typing
               └─> calculate_multiomics_sample_embedding
                     └─> (optional) find_optimal_cell_resolution_multiomics_linux
                            └─> replace_optimal_dimension_reduction
                                  └─> downstream analyses
```

For a runnable walkthrough, see the [Multi-omics tutorial](../tutorials/multiomics.md).
