# Downstream API

This page documents downstream functions called after sample embeddings are available.

## `sample_distance(...)`

Source: `code/sample_distance/sample_distance.py`

```python
sample_distance(
    adata,
    output_dir,
    method,
    data_type="ATAC",
    grouping_columns=None,
    summary_csv_path=None,
    cell_adata=None,
    cell_type_column="cell_type",
    sample_column="sample",
    embedding_key="X_pca_harmony",
    n_pcs=20,
    proportions=None,
    centroids=None,
    pseudobulk_adata=None,
)
```

Compute sample distance outputs (cosine/correlation/EMD/Jensen-Shannon/Chi-square).

## `CCA_Call(...)`

Source: `code/sample_trajectory/CCA.py`

```python
CCA_Call(
    adata,
    output_dir=None,
    trajectory_col="sev.level",
    n_components=2,
    auto_select_best_2pc=True,
    verbose=False,
    show_sample_labels=False,
)
```

Run supervised CCA-based trajectory inference.

## `cca_pvalue_test(...)`

Source: `code/sample_trajectory/CCA_test.py`

```python
cca_pvalue_test(
    pseudo_adata,
    column,
    input_correlation,
    output_directory,
    num_simulations=1000,
    trajectory_col="sev.level",
    verbose=True,
)
```

Compute empirical p-value diagnostics for CCA trajectory correlation.

## `TSCAN(...)`

Source: `code/sample_trajectory/TSCAN.py`

```python
TSCAN(
    AnnData_sample,
    column,
    n_clusters=None,
    output_dir="./",
    grouping_columns=None,
    verbose=False,
    origin=None,
    pseudotime_mode="rank",
)
```

Run unsupervised sample-level pseudotime inference.

## `run_trajectory_gam_differential_gene_analysis(...)`

Source: `code/sample_trajectory/trajectory_diff_gene.py`

```python
run_trajectory_gam_differential_gene_analysis(
    pseudobulk_adata,
    pseudotime_source,
    *,
    sample_col="sample",
    pseudotime_col="pseudotime",
    covariate_columns=None,
    fdr_threshold=0.01,
    effect_size_threshold=1.0,
    top_n_genes=100,
    num_splines=5,
    spline_order=3,
    output_dir="trajectory_diff_gene_results_single",
    visualization_gene_list=None,
    generate_visualizations=True,
    group_col=None,
    n_clusters=3,
    top_n_genes_for_curves=20,
    verbose=True,
)
```

Identify trajectory-associated genes/features using GAM models.

## `cluster(...)`

Source: `code/cluster.py`

```python
cluster(
    pseudobulk_adata,
    output_dir,
    number_of_clusters=5,
    use_expression=True,
    use_proportion=True,
    random_state=0,
)
```

K-means clustering in expression/proportion sample embedding spaces.

## `proportion_test(...)`

Source: `code/sample_clustering/proportion_test.py`

```python
proportion_test(
    adata,
    sample_col,
    group_col=None,
    sample_to_clade=None,
    celltype_col="celltype",
    output_dir=None,
    verbose=True,
)
```

Test cell type composition differences across groups or clusters.

## `raisinfit(...)`

Source: `code/sample_clustering/RAISIN.py`

```python
raisinfit(
    adata,
    sample_col,
    testtype="unpaired",
    group_col=None,
    individual_col=None,
    batch_col=None,
    sample_to_clade=None,
    custom_design=None,
    intercept=True,
    filtergene=False,
    filtergenequantile=0.5,
    n_jobs=None,
    verbose=True,
)
```

Fit RAISIN model for cluster differential analysis.

## `run_pairwise_tests(...)`

Source: `code/sample_clustering/RAISIN_TEST.py`

```python
run_pairwise_tests(
    fit,
    output_dir,
    groups_to_compare=None,
    control_group=None,
    fdrmethod="fdr_bh",
    n_permutations=100,
    fdr_threshold=0.05,
    top_n_genes=50,
    make_summary_plots=True,
    verbose=True,
)
```

Run pairwise RAISIN contrasts and produce summary plots/tables.

## `visualization(...)`

Source: `code/visualization/visualization_other.py`

```python
visualization(
    AnnData_cell,
    pseudobulk_anndata,
    output_dir,
    grouping_columns=None,
    age_bin_size=None,
    age_column="age",
    verbose=True,
    plot_dendrogram_flag=True,
    plot_cell_type_proportions_pca_flag=False,
    plot_cell_type_expression_umap_flag=False,
)
```

Generate downstream visualization panels (dendrogram and optional embeddings).

## `visualize_multimodal_embedding(...)`

Source: `code/visualization/multi_omics_visualization.py`

```python
visualize_multimodal_embedding(
    adata,
    modality_col=None,
    color_col=None,
    target_modality=None,
    expression_key="X_DR_expression",
    proportion_key="X_DR_proportion",
    figsize=(20, 8),
    point_size=60,
    alpha=0.8,
    colormap="viridis",
    categorical_cmap="tab10",
    output_dir=None,
    show_sample_names=False,
    force_data_type=None,
    show_default=True,
    verbose=True,
    visualization_grouping_column=None,
)
```

Visualize multi-omics sample embeddings with modality/group coloring.
