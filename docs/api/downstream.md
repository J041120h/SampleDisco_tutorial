# API Reference: Downstream Analysis

The shared `downstream_analysis(...)` function in `code/wrapper/wrapper.py` dispatches the major sample-level analyses after embeddings are computed.

## `downstream_analysis(...)`

```python
downstream_analysis(
    pseudo_adata,
    output_dir: str,
    modality: str,
    status_flags: dict,
    sample_distance_calculation: bool = True,
    trajectory_analysis: bool = True,
    trajectory_DGE: bool = True,
    sample_cluster: bool = True,
    proportion_test: bool = False,
    cluster_DGE: bool = False,
    visualize_data: bool = True,
    ...
) -> dict
```

**Responsibilities**

- Run sample distance calculations.
- Perform supervised CCA or unsupervised TSCAN trajectory inference.
- Launch GAM-based trajectory differential testing.
- Cluster samples and run optional group comparisons.
- Produce auxiliary visualizations.

## `sample_distance(...)`

```python
sample_distance(
    adata: AnnData,
    output_dir: str,
    method: str,
    data_type: str = "ATAC",
    ...
)
```

Computes sample-sample distances from reduced sample embeddings or cell-type composition summaries.

**Supported families**

- Standard distances such as cosine and correlation on reduced sample coordinates.
- EMD using cell-type proportions and centroids.
- Chi-square and Jensen-Shannon distances on compositional profiles.

## `CCA_Call(...)` and `TSCAN(...)`

```python
CCA_Call(adata: AnnData, output_dir: str = None, trajectory_col: str = "sev.level", ...)
TSCAN(AnnData_sample: sc.AnnData, column: str, output_dir: str = "./", ...)
```

Two trajectory modes are exposed:

- `CCA_Call(...)` for supervised trajectory alignment when a known ordering column is available.
- `TSCAN(...)` for unsupervised pseudotime estimation from sample embeddings.

## `run_trajectory_gam_differential_gene_analysis(...)`

```python
run_trajectory_gam_differential_gene_analysis(
    pseudobulk_adata: ad.AnnData,
    pseudotime_source,
    sample_col: str = "sample",
    pseudotime_col: str = "pseudotime",
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 1,
    top_n_genes: int = 100,
    ...
) -> pd.DataFrame
```

Fits GAM-style models along pseudotime and returns a table of significant trajectory-associated genes or features.

## `cluster(...)`

```python
cluster(
    pseudobulk_adata: ad.AnnData,
    output_dir: str,
    number_of_clusters: int = 5,
    use_expression: bool = True,
    use_proportion: bool = True,
    random_state: int = 0,
) -> Tuple[Optional[Dict[str, int]], Optional[Dict[str, int]]]
```

Clusters samples in expression-derived and/or proportion-derived spaces.

## `proportion_test(...)`, `raisinfit(...)`, and `run_pairwise_tests(...)`

These functions perform group-aware statistical testing on sample clusters and cell-type proportions.

- `proportion_test(...)`: assess composition differences.
- `raisinfit(...)`: fit the RAISIN model for paired or unpaired comparisons.
- `run_pairwise_tests(...)`: extract pairwise group contrasts and result tables.

## Visualization helpers

- `visualization(...)` in `visualization/visualization_other.py`
- `visualize_multimodal_embedding(...)` in `visualization/multi_omics_visualization.py`

Use these after sample embeddings are available to generate PCA, UMAP, dendrogram, or modality-colored summaries.

## Suggested reading order

1. [Preparation API](preparation.md)
2. [Sample Embedding API](embedding.md)
3. [Tutorial pages](../tutorials/config_overview.md)
