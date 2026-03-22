# Downstream API

## `downstream_analysis(...)`

```python
downstream_analysis(
    pseudo_adata,
    output_dir: str,
    modality: str,
    status_flags: dict,
    adata_cell=None,
    adata_sample=None,
    sample_distance_calculation: bool = True,
    trajectory_analysis: bool = True,
    trajectory_DGE: bool = True,
    sample_cluster: bool = True,
    proportion_test: bool = False,
    cluster_DGE: bool = False,
    visualize_data: bool = True,
    visualize_embedding: bool = False,
    use_gpu: bool = False,
    verbose: bool = True,
    sample_col: str = "sample",
    batch_col: Optional[Union[str, List[str]]] = None,
    celltype_col: str = "cell_type",
    sample_distance_methods: Optional[List[str]] = None,
    grouping_columns: Optional[List[str]] = None,
    summary_sample_csv_path: Optional[str] = None,
    n_cca_pcs: int = 2,
    trajectory_col: str = "sev.level",
    trajectory_supervised: bool = False,
    trajectory_visualization_label: Optional[List[str]] = None,
    cca_pvalue: bool = False,
    tscan_origin: Optional[str] = None,
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 1,
    top_n_genes: int = 100,
    trajectory_diff_gene_covariate: Optional[List] = None,
    num_splines: int = 5,
    spline_order: int = 3,
    visualization_gene_list: Optional[List] = None,
    cluster_number: int = 4,
    cluster_differential_gene_group_col: Optional[str] = None,
    age_bin_size: Optional[int] = None,
    age_column: str = "age",
    plot_dendrogram_flag: bool = True,
    plot_cell_type_proportions_pca_flag: bool = False,
    plot_cell_type_expression_umap_flag: bool = False,
    multiomics_modality_col: str = "modality",
    multiomics_color_col: Optional[str] = None,
    multiomics_visualization_grouping_column: Optional[List[str]] = None,
    multiomics_target_modality: str = "ATAC",
    multiomics_expression_key: str = "X_DR_expression",
    multiomics_proportion_key: str = "X_DR_proportion",
    multiomics_figsize: Tuple[int, int] = (20, 8),
    multiomics_point_size: int = 60,
    multiomics_alpha: float = 0.8,
    multiomics_colormap: str = "viridis",
    multiomics_show_sample_names: bool = False,
    multiomics_force_data_type: Optional[str] = None,
) -> dict
```

Shared orchestrator for downstream modules across RNA/ATAC/multi-omics.

### Main controls

| Group | Key parameters | Meaning |
| --- | --- | --- |
| Distance | `sample_distance_calculation`, `sample_distance_methods`, `grouping_columns` | Sample-to-sample distance generation |
| Trajectory | `trajectory_analysis`, `trajectory_supervised`, `trajectory_col`, `n_cca_pcs` | CCA/TSCAN pseudotime inference |
| Trajectory DGE | `trajectory_DGE`, `fdr_threshold`, `effect_size_threshold`, `num_splines` | GAM-based trajectory differential analysis |
| Clustering | `sample_cluster`, `cluster_number` | K-means clustering on sample embeddings |
| Proportion/cluster DGE | `proportion_test`, `cluster_DGE`, `cluster_differential_gene_group_col` | Proportion testing and RAISIN analysis |
| Visualization | `visualize_data`, `plot_*` | Dendrogram and embedding visualization controls |
| Multi-omics embedding viz | `visualize_embedding`, `multiomics_*` | Multi-modal embedding display controls |

## `sample_distance(...)`

- Source: `sample_distance/sample_distance.py`
- Purpose: compute sample distances from expression/proportion DR or compositional methods.
- Important `method` values include `cosine`, `correlation`, `emd`, `jensenshannon`, and `chi_square`.

## `CCA_Call(...)`, `cca_pvalue_test(...)`, `TSCAN(...)`

- `CCA_Call()` (source: `sample_trajectory/CCA.py`) runs supervised CCA trajectory.
- `cca_pvalue_test()` (source: `sample_trajectory/CCA_test.py`) runs p-value checks for CCA correlation.
- `TSCAN()` (source: `sample_trajectory/TSCAN.py`) runs unsupervised pseudotime.

## `run_trajectory_gam_differential_gene_analysis(...)`

- Source: `sample_trajectory/trajectory_diff_gene.py`
- Purpose: fit GAM models across pseudotime and detect trajectory-associated genes/features.

## `cluster(...)`

- Source: `cluster.py`
- Purpose: K-means clustering for expression/proportion sample embeddings.

## `proportion_test(...)`

- Source: `sample_clustering/proportion_test.py`
- Purpose: statistical testing on cell type composition differences across groups/clusters.

## `raisinfit(...)` and `run_pairwise_tests(...)`

- Sources: `sample_clustering/RAISIN.py`, `sample_clustering/RAISIN_TEST.py`
- Purpose: cluster differential analysis and pairwise tests.

## `visualization(...)` and `visualize_multimodal_embedding(...)`

- `visualization()` source: `visualization/visualization_other.py`
- `visualize_multimodal_embedding()` source: `visualization/multi_omics_visualization.py`

Use these to generate standard downstream visual outputs after embeddings exist.

## Related

- [RNA API](rna.md)
- [ATAC API](atac.md)
- [Multi-omics API](multiomics.md)
- [Downstream tutorial](../tutorials/tutorial_downstream.md)
