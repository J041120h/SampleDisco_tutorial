# Downstream analysis

Once the sample embedding has been computed, every downstream module runs off the same object. The sample embedding is a single key, `adata.uns['X_DR_sample']` (a `samples × PCs` DataFrame), written in place by `compute_sample_embedding` into the cell-level AnnData. The functions are identical across RNA, ATAC, and multi-omics — only the input AnnData changes.

!!! note "Figures on these pages use the RNA COVID run"
    All example figures shown across the downstream tutorials come from the scRNA-seq COVID pipeline output. The calls are the same for ATAC and multi-omics; point the functions at the matching AnnData and the downstream behavior is identical.

## Before you start

The cell-level AnnData written by preprocessing (`adata_preprocessed.h5ad`) carries the sample embedding in `uns['X_DR_sample']` after `compute_sample_embedding` runs. Most downstream modules take a one-row-per-sample AnnData; build it on the fly from that key with `build_sample_adata`.

```python
import anndata as ad
from sampledisco.sample_embedding import compute_sample_embedding
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

adata_cell = ad.read_h5ad("/results/rna/preprocess/adata_preprocessed.h5ad")

# Compute the sample embedding in place (skip if uns['X_DR_sample'] already present)
adata_cell  = compute_sample_embedding(adata_cell, output_dir="/results/rna")

# Materialize the sample-level AnnData (samples × PCs) for downstream consumers
pseudo_adata = build_sample_adata(adata_cell, sample_col="sample")
```

!!! tip "The config-driven wrapper does all of this for you"
    Running the full pipeline with `sampledisco -m complex --config <config.yaml>` preprocesses, computes the sample embedding, and runs every downstream module below in one call — you rarely need to invoke these functions by hand.

## Tutorials

Each function has its own page with the call, parameters, outputs, and representative figures.

| Page | Function | What it does |
| --- | --- | --- |
| [Resolution selection](resolution_selection.md) | _(deprecated)_ | The CCA-guided resolution sweep was removed; superseded by alpha / block-weight autotune (`run_autotune`, the `*_autotune_enable` wrapper flags). |
| [Sample distance](sample_distance.md) | `sample_distance` | Pairwise distance matrices on the sample embedding with multiple metrics. |
| [Dimension association](sample_association.md) | `run_dimension_association_analysis` | Per-PC variance-explained decomposition against every metadata variable — confounder / leading-covariate check. |
| [Trajectory — CCA](trajectory_cca.md) | `CCA_Call`, `cca_pvalue_test` | Supervised pseudotime and its permutation significance. |
| [Trajectory — TSCAN](trajectory_tscan.md) | `TSCAN` | Unsupervised pseudotime via GMM + MST. |
| [Trajectory DGE](trajectory_dge.md) | `run_trajectory_gam_differential_gene_analysis` | GAM-based differential expression along pseudotime. |
| [Sample clustering](cluster.md) | `cluster` | K-means on the sample embedding. |
| [Proportion test](proportion_test.md) | `proportion_test` | Cell-type composition test across sample groups. |
| [RAISIN cluster DGE](raisin_dge.md) | `raisinfit`, `run_pairwise_tests` | Hierarchical GLM for cluster-level differential expression. |
| [General visualization](visualization.md) | `visualization` | Cell-type dendrogram, proportion PCA, expression UMAP. |

Recommended reading order: sample distance → dimension association → trajectory (CCA or TSCAN) → trajectory DGE → sample clustering → proportion test / RAISIN → general visualization.
