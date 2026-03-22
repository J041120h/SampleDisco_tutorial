# Downstream Analysis Tutorial

This page documents the shared `downstream_analysis(...)` stage used by RNA, ATAC, and multi-omics wrappers. Example outputs below use RNA results.

## Main function call

```python
from wrapper.wrapper import downstream_analysis

downstream_analysis(
    pseudo_adata=pseudo_adata,
    output_dir=".../rna",
    modality="rna",
    status_flags=status_flags,
    sample_distance_calculation=True,
    trajectory_analysis=True,
    trajectory_DGE=True,
    sample_cluster=True,
    proportion_test=True,
    cluster_DGE=False,
    visualize_data=True,
)
```

## 1) Sample distance

Main API: `sample_distance()`

![Distance expression cosine](../assets/images/downstream/sample_distance_expression_DR_heatmap_cosine.pdf)
![Distance proportion cosine](../assets/images/downstream/sample_distance_proportion_DR_heatmap_cosine.pdf)
![Distance expression correlation](../assets/images/downstream/sample_distance_expression_DR_heatmap_correlation.pdf)
![Distance proportion correlation](../assets/images/downstream/sample_distance_proportion_DR_heatmap_correlation.pdf)

## 2) Supervised trajectory (CCA)

Main API: `CCA_Call()`, optional `cca_pvalue_test()`

![CCA expression](../assets/images/downstream/cca_expression.pdf)
![CCA expression contributions](../assets/images/downstream/cca_expression_contributions.pdf)
![CCA proportion](../assets/images/downstream/cca_proportion.pdf)
![CCA p-value expression](../assets/images/downstream/cca_pvalue_X_DR_expression.png)
![CCA p-value proportion](../assets/images/downstream/cca_pvalue_X_DR_proportion.png)

## 3) Unsupervised trajectory (TSCAN)

Main API: `TSCAN()`

![TSCAN clusters by cluster expression](../assets/images/downstream/tscan_clusters_by_cluster_expression.png)
![TSCAN clusters by grouping expression](../assets/images/downstream/tscan_clusters_by_grouping_expression.png)
![TSCAN clusters by cluster proportion](../assets/images/downstream/tscan_clusters_by_cluster_proportion.png)
![TSCAN clusters by grouping proportion](../assets/images/downstream/tscan_clusters_by_grouping_proportion.png)

## 4) Trajectory DGE (GAM)

Main API: `run_trajectory_gam_differential_gene_analysis()`

![Trajectory DGE summary](../assets/images/downstream/results_summary.png)
![Trajectory DGE heatmap](../assets/images/downstream/tde_heatmap.png)
![Trajectory DGE volcano](../assets/images/downstream/volcano_plot.png)
![Trajectory DGE gene curves](../assets/images/downstream/gene_curves.png)
![Trajectory DGE sample density](../assets/images/downstream/sample_density.png)
![Trajectory DGE cluster patterns](../assets/images/downstream/cluster_patterns.png)
![Trajectory DGE sample-level curves](../assets/images/downstream/sample_level_curves.png)

## 5) Sample clustering

Main API: `cluster()`

![Kmeans expression embedding](../assets/images/downstream/kmeans_expression_embedding.png)
![Kmeans proportion embedding](../assets/images/downstream/kmeans_proportion_embedding.png)

## 6) Proportion test

Main API: `proportion_test()`

![Proportion heatmap](../assets/images/downstream/proportion_heatmap_group_by_celltype.png)
![Proportion significance matrix](../assets/images/downstream/proportion_significance_matrix.png)

## 7) Cluster DGE (RAISIN)

Main API: `raisinfit()` and `run_pairwise_tests()`

Enable this module with:

- `*_cluster_dge: true`
- `*_cluster_differential_gene_group_col`

## 8) Additional visualization

Main API: `visualization()`

Use these controls:

- `*_plot_dendrogram_flag`
- `*_plot_cell_type_proportions_pca_flag`
- `*_plot_cell_type_expression_umap_flag`

