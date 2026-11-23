# `downstream_analysis`

Phase 2 of the [`wrapper`](wrapper.md): runs the downstream steps on a finished sample embedding, with fixed output folders. Call it directly to re-run or add downstream steps on an existing embedding without repeating preprocessing, cell typing and the embedding. Steps run in this order, each behind a flag:

| Step (flag) | Function | Output under `output_dir` |
| --- | --- | --- |
| `sample_distance_calculation` | [`sample_distance`](../downstream/sample_distance.md), once per method in `sample_distance_methods` | `Sample_distance/<method>/` |
| `trajectory_analysis` | `trajectory_supervised=True`: [`CCA_Call`](../downstream/trajectory_cca_call.md) on `n_cca_pcs` PCs, plus [`cca_pvalue_test`](../downstream/trajectory_cca_pvalue_test.md) with `cca_pvalue=True`. Otherwise [`TSCAN`](../downstream/trajectory_tscan.md). The pseudotime is stored in `pseudo_adata.obs['pseudotime_sample']`. | `CCA/`, `CCA_test/` or `TSCAN/` |
| `trajectory_DGE` (only inside `trajectory_analysis`) | [`run_trajectory_gam_differential_gene_analysis`](../downstream/trajectory_gam_dge.md) on `adata_cell_for_dge` or `adata_cell`; with TSCAN, `anchor_col=trajectory_col` | `trajectoryDEG/` |
| `sample_cluster` | [`cluster`](../downstream/cluster.md) with `cluster_number` (seed 0) | `sample_cluster/` |
| `proportion_test` | [`proportion_test`](../downstream/proportion_test.md) on `adata_sample`, grouped by `cluster_differential_gene_group_col` if set, else by the k-means clusters | `sample_cluster/proportion_test/` |
| `cluster_DGE` | [`raisinfit`](../downstream/raisinfit.md) (unpaired, first `batch_col`, same grouping) + [`run_pairwise_tests`](../downstream/run_pairwise_tests.md) | `raisin_results/` |
| `visualize_data` | [`visualization`](../downstream/visualization.md) | `visualization/` |
| `visualize_embedding` (`modality="multiomics"` only) | [`visualize_multimodal_embedding`](../downstream/visualize_multimodal_embedding.md) in its default mode | `visualization/` |
| `dimension_association_analysis` | [`run_dimension_association_analysis`](../downstream/sample_association.md) | `sample_association/` |

Errors in `proportion_test`, `cluster_DGE`, `visualize_embedding` and `dimension_association_analysis` are printed and skipped; errors in the other steps are raised.

**Import:** `from sampledisco.wrapper.wrapper import downstream_analysis`

**Source:** `wrapper/wrapper.py:100`

## Signature

```python
def downstream_analysis(
    pseudo_adata,
    output_dir: str,
    modality: str,
    status_flags: dict,
    adata_cell=None,
    adata_cell_for_dge=None,
    adata_sample_for_dge=None,
    adata_sample=None,
    dge_modality: Optional[str] = None,
    # step flags
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
    random_state: int = 42,
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
    tscan_n_clusters: Optional[int] = None,
    tscan_pseudotime_mode: str = "rank",
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
    dimension_association_analysis: bool = False,
    association_continuous_cols: Optional[List[str]] = None,
    association_categorical_cols: Optional[List[str]] = None,
    association_n_permutations: int = 999,
    multiomics_modality_col: str = "modality",
    multiomics_color_col: Optional[str] = None,
    multiomics_visualization_grouping_column: Optional[List[str]] = None,
    multiomics_target_modality: str = "ATAC",
    multiomics_sample_embedding_key: str = "X_DR_sample",
    multiomics_figsize: Tuple[int, int] = (20, 8),
    multiomics_point_size: int = 60,
    multiomics_alpha: float = 0.8,
    multiomics_colormap: str = "viridis",
    multiomics_show_sample_names: bool = False,
    multiomics_force_data_type: Optional[str] = None,
    dge_pseudobulk_celltype_col: Optional[str] = None,
    dge_pseudobulk_batch_col: Optional[Union[str, List[str]]] = None,
    dge_pseudobulk_n_features_per_celltype: Optional[int] = 2000,
    dge_pseudobulk_columns_to_preserve: Optional[Union[str, List[str]]] = None,
) -> dict
```

## Parameters

The step parameters have the same meaning as the wrapper's `rna_*` / `atac_*` / `multiomics_*` keys without the prefix (`rna_cluster_number` → `cluster_number`, …); see the [config key reference](../../tutorials/config_reference.md) and the linked function pages. The ones specific to a direct call:

| Name | Default | Description |
| --- | --- | --- |
| `pseudo_adata` | — | Sample-level AnnData from [`build_sample_adata`](../shared/build_sample_adata.md). `None` skips every step and returns at once. |
| `output_dir` | — | Root folder for the outputs in the table above. |
| `modality` | — | `"rna"`, `"atac"` or `"multiomics"`. Sets the distance `data_type` and enables `visualize_embedding` for `"multiomics"`. |
| `status_flags` | — | A dict that **must contain the key `modality`**, e.g. `{"rna": {}}`. Completed steps are recorded in it. The proportion / expression scatter plots of `visualize_data` additionally require `status_flags[modality]["derive_sample_embedding"]` to be true. |
| `adata_cell` | `None` | Cell-level AnnData: cell types for the dendrogram (required when `visualize_data` and `plot_dendrogram_flag` are on), cells for EMD / chi-square / Jensen–Shannon distances, and DGE input. |
| `adata_sample` | `None` | Cell-level AnnData for `proportion_test` and RAISIN (the wrapper passes the same object as `adata_cell`). |
| `adata_cell_for_dge` / `adata_sample_for_dge` / `dge_modality` | `None` | Multi-omics: the RNA-only file with real expression for DGE / RAISIN; `dge_modality="RNA"` maps `<sample>_RNA` units back to bare sample ids. |
| `grouping_columns` | `None` → `["sev.level"]` | Grouping for the distance check and the plots. |
| `sample_distance_methods` | `None` → `["cosine", "correlation"]` | Distance methods. |
| `trajectory_visualization_label` | `None` → `["sev.level"]` | TSCAN plot overlays. |
| `summary_sample_csv_path` | `None` → `{output_dir}/summary_sample.csv` | Group-separation summary filled by the distance checks. A file created by this call is **deleted** after the trajectory step. |
| `use_gpu` | `False` | Accepted; not used by any downstream step. |
| `random_state` | `42` | Seed for RAISIN. |

## Returns

`{"pseudo_adata": pseudo_adata, "status_flags": status_flags}`. `pseudo_adata` is the input object, now with `obs['pseudotime_sample']` (after a trajectory step) and `obs['cluster_sample_kmeans']` (after clustering).

## Usage

Continuing from the RNA chain: sample distances and k-means only.

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.wrapper.wrapper import downstream_analysis

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

out = downstream_analysis(
    pseudo,
    "api_demo/rna_downstream",
    modality="rna",
    status_flags={"rna": {}},
    adata_cell=adata,
    adata_sample=adata,
    sample_distance_methods=["cosine"],
    grouping_columns=["sev.level"],
    trajectory_analysis=False,
    visualize_data=False,
    cluster_number=2,
)
print(out["status_flags"])
print(out["pseudo_adata"].obs["cluster_sample_kmeans"].to_dict())
```

**Writes** → under `api_demo/rna_downstream/`: `Sample_distance/cosine/` (the five files listed on [`sample_distance`](../downstream/sample_distance.md)) and `sample_cluster/kmeans_clusters_sample.csv`, `kmeans_sample_embedding.png`. The temporary `summary_sample.csv` is created and deleted again. Printed (about 2 s):

```text
{'rna': {'sample_distance_calculation': True, 'sample_cluster': True}}
{'CoV-3-Wilk': '0', 'CoV-34-Aruna': '1', 'CoV-36-Aruna': '0', 'CoV-4-Wilk': '1', 'HD-1-Wilk': '0', 'HD-20-Aruna': '0', 'HD-21-Aruna': '0', 'HD-4-Wilk': '0'}
```
