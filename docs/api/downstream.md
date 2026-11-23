# Downstream API

Everything in this section runs on the sample embedding `X_DR_sample` (units × PCs) written by [`compute_sample_embedding`](shared/calculate_sample_embedding.md) or [`run_autotune`](shared/run_autotune.md). Most functions take the **sample-level** AnnData from [`build_sample_adata`](shared/build_sample_adata.md) (`pseudo` in the examples); `proportion_test`, `raisinfit` and the trajectory DGE take the **cell-level** AnnData (`adata`). The functions are the same for RNA, ATAC and multi-omics. In a pipeline run, [`downstream_analysis`](pipeline/downstream_analysis.md) calls them with fixed output folders.

The examples on these pages start from the RNA chain's output; see [Running the examples](index.md#running-the-examples).

## Sample distance

| Function | Purpose |
| --- | --- |
| [`sample_distance`](downstream/sample_distance.md) | Pairwise distances between samples: any `scipy` `pdist` metric on the embedding (cosine, correlation, euclidean, …), or EMD, chi-square and Jensen–Shannon on cell-type proportions. Optionally scores between- vs within-group distance. |

## Dimension association

| Function | Purpose |
| --- | --- |
| [`run_dimension_association_analysis`](downstream/sample_association.md) | R² of each embedding PC against each sample-level variable, with permutation p-values: which covariates drive which PCs. |

## Trajectory analysis

| Function | Purpose | When to use |
| --- | --- | --- |
| [`CCA_Call`](downstream/trajectory_cca_call.md) | Supervised pseudotime: CCA between the leading PCs and a numeric phenotype; samples are ordered along the CCA direction. | You have an ordered phenotype (severity, age, stage). |
| [`cca_pvalue_test`](downstream/trajectory_cca_pvalue_test.md) | Permutation p-value for a CCA correlation. | Testing the CCA result. |
| [`TSCAN`](downstream/trajectory_tscan.md) | Unsupervised pseudotime: Gaussian-mixture clusters, minimum spanning tree, projection on the main path. | No phenotype to supervise with. |
| [`run_trajectory_gam_differential_gene_analysis`](downstream/trajectory_gam_dge.md) | Per-gene GAM of pseudobulk expression against pseudotime. | Genes that change along a CCA or TSCAN pseudotime. |

## Sample clustering and group comparisons

| Function | Purpose | When to use |
| --- | --- | --- |
| [`cluster`](downstream/cluster.md) | K-means on the embedding. | A discrete grouping of samples. |
| [`proportion_test`](downstream/proportion_test.md) | CLR + limma-style moderated test of cell-type proportions between groups. | Is a group difference driven by composition? |
| [`raisinfit`](downstream/raisinfit.md) | Fit the RAISIN hierarchical model on cell-level expression. | Step 1 of group-level differential expression. |
| [`run_pairwise_tests`](downstream/run_pairwise_tests.md) | All pairwise RAISIN contrasts, with volcano and summary plots. | Step 2 of group-level differential expression. |

## Prediction and label transfer

| Function | Purpose | When to use |
| --- | --- | --- |
| [`predict_phenotype`](downstream/predict_phenotype.md) | Cross-validated logistic or ridge regression on the PCs, with donor-grouped folds. | How well does the embedding predict a phenotype? |
| [`transfer_labels`](downstream/transfer_labels.md) | k-nearest-neighbor label transfer with a leave-one-out accuracy. | Label new samples from labelled ones. |

## Visualization

| Function | Purpose |
| --- | --- |
| [`visualization`](downstream/visualization.md) | Cell-type dendrogram (Ward linkage on mean `Z_comp` per cell type) and scatter plots of the first two embedding PCs colored by each grouping column. |
| [`visualize_multimodal_embedding`](downstream/visualize_multimodal_embedding.md) | Scatter plots of the first two embedding PCs for multi-omics units, colored by columns or with one modality highlighted. |

!!! note "Eight samples"
    On the 8-sample demo, CCA correlations, permutation p-values, prediction scores and DGE calls are degenerate or noisy; the examples show the calls and the files they write. The limits are listed under [Small cohorts](../tutorials/downstream_details.md#small-cohorts).

For a worked run of every analysis, see [Downstream analysis](../tutorials/downstream/index.md); for options and caveats, [Downstream options and caveats](../tutorials/downstream_details.md).
