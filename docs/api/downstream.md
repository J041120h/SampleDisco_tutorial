# Downstream API

Everything in this section runs off a `pseudobulk_sample.h5ad` — the sample-level AnnData with `X_DR_expression` and `X_DR_proportion` in `.obsm`. The functions are identical across RNA, ATAC, and multi-omics, and most share a common wrapper invocation inside `downstream_analysis()`.

## Sample distance

| Function | Purpose |
| --- | --- |
| [`sample_distance`](downstream/sample_distance.md) | Unified entry point for vector metrics (cosine, correlation, euclidean, ...), EMD, chi-square, Jensen–Shannon. |

## Trajectory analysis

| Function | Purpose | When to use |
| --- | --- | --- |
| [`CCA_Call`](downstream/trajectory_cca_call.md) | Supervised pseudotime via canonical correlation with a phenotype. | You have a phenotype gradient in `.obs`. |
| [`cca_pvalue_test`](downstream/trajectory_cca_pvalue_test.md) | Permutation p-value for a CCA correlation. | Want to test significance of a CCA result. |
| [`TSCAN`](downstream/trajectory_tscan.md) | Unsupervised pseudotime via GMM + MST. | No phenotype; let structure emerge from the embedding. |
| [`run_trajectory_gam_differential_gene_analysis`](downstream/trajectory_gam_dge.md) | Per-gene GAM against pseudotime, with heatmap, volcano, and curve plots. | Pseudotime (from CCA or TSCAN) in hand; find trajectory-associated genes. |

## Sample clustering and cluster comparisons

| Function | Purpose | When to use |
| --- | --- | --- |
| [`cluster`](downstream/cluster.md) | K-means on the two sample embeddings. | Want a discrete sample grouping when no phenotype is available. |
| [`proportion_test`](downstream/proportion_test.md) | eBayes moderated test on per-sample cell-type proportions across groups. | Ask whether cluster/group differences are driven by cell-type composition. |
| [`raisinfit`](downstream/raisinfit.md) | Fit the RAISIN hierarchical GLM. | Step 1 of cluster-level differential expression. |
| [`run_pairwise_tests`](downstream/run_pairwise_tests.md) | Run every pairwise contrast on a fitted RAISIN model, emit volcanos + summary plots. | Step 2 of cluster-level DE. |

## Visualization

| Function | Purpose |
| --- | --- |
| [`visualization`](downstream/visualization.md) | Cell-type dendrogram, proportion PCA, expression UMAP — toggle which panels you want. |
| [`visualize_multimodal_embedding`](downstream/visualize_multimodal_embedding.md) | Side-by-side scatter of the two embeddings with flexible coloring; multi-omics only. |

For an end-to-end walkthrough that chains these functions together, see the [Downstream tutorials](../tutorials/downstream/index.md).
