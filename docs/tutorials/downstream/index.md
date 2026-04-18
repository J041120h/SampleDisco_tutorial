# Downstream analysis

Once you have a `pseudobulk_sample.h5ad` from any of the three pipelines, every downstream module runs off the same object. The functions are identical across RNA, ATAC, and multi-omics — only the input AnnData changes.

!!! note "Figures on these pages use the RNA COVID run"
    All example figures shown across the downstream tutorials come from the scRNA-seq COVID pipeline output. The calls are the same for ATAC and multi-omics; point the functions at the matching `pseudo_adata` and the downstream behavior is identical.

## Before you start

```python
import anndata as ad

pseudo_adata = ad.read_h5ad("/results/rna/pseudobulk/pseudobulk_sample.h5ad")
adata_cell   = ad.read_h5ad("/results/rna/preprocess/adata_cell.h5ad")
```

## Tutorials

Each function has its own page with the call, parameters, outputs, and representative figures.

| Page | Function | What it does |
| --- | --- | --- |
| [Resolution selection](resolution_selection.md) | `find_optimal_cell_resolution_linux` | CCA-guided two-pass search for the Leiden resolution that best aligns with a phenotype. |
| [Sample distance](sample_distance.md) | `sample_distance` | Pairwise distance matrices on either embedding with multiple metrics. |
| [Trajectory — CCA](trajectory_cca.md) | `CCA_Call`, `cca_pvalue_test` | Supervised pseudotime and its permutation significance. |
| [Trajectory — TSCAN](trajectory_tscan.md) | `TSCAN` | Unsupervised pseudotime via GMM + MST. |
| [Trajectory DGE](trajectory_dge.md) | `run_trajectory_gam_differential_gene_analysis` | GAM-based differential expression along pseudotime. |
| [Sample clustering](cluster.md) | `cluster` | K-means on each sample embedding. |
| [Proportion test](proportion_test.md) | `proportion_test` | Cell-type composition test across sample groups. |
| [RAISIN cluster DGE](raisin_dge.md) | `raisinfit`, `run_pairwise_tests` | Hierarchical GLM for cluster-level differential expression. |
| [General visualization](visualization.md) | `visualization` | Cell-type dendrogram, proportion PCA, expression UMAP. |

Recommended reading order: resolution selection (optional) → sample distance → trajectory (CCA or TSCAN) → trajectory DGE → sample clustering → proportion test / RAISIN → general visualization.
