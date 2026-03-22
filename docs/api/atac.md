# ATAC API

## `atac_wrapper(...)`

```python
atac_wrapper(
    atac_count_data_path: str = None,
    atac_output_dir: str = None,
    preprocessing: bool = True,
    cell_type_cluster: bool = True,
    derive_sample_embedding: bool = True,
    cca_based_cell_resolution_selection: bool = False,
    use_gpu: bool = False,
    verbose: bool = True,
    status_flags: dict = None,
    adata_cell_path: str = None,
    adata_sample_path: str = None,
    atac_sample_meta_path: str = None,
    cell_meta_path: str = None,
    pseudo_adata_path: str = None,
    sample_col: str = "sample",
    sample_level_batch_col: Optional[List[str]] = None,
    celltype_col: str = "cell_type",
    cell_embedding_column: str = None,
    min_cells: int = 1,
    min_features: int = 2000,
    max_features: int = 15000,
    min_cells_per_sample: int = 1,
    exclude_features: list = None,
    cell_level_batch_key: list = None,
    doublet_detection: bool = True,
    num_cell_hvfs: int = 50000,
    cell_embedding_num_pcs: int = 50,
    num_harmony_iterations: int = 30,
    tfidf_scale_factor: float = 1e4,
    log_transform: bool = True,
    drop_first_lsi: bool = True,
    leiden_cluster_resolution: float = 0.8,
    existing_cell_types: bool = False,
    n_target_cell_clusters: int = None,
    umap: bool = False,
    sample_hvg_number: int = 50000,
    sample_embedding_dimension: int = 30,
    harmony_for_proportion: bool = True,
    preserve_cols_in_sample_embedding: list = None,
    trajectory_col: str = "sev.level",
    n_cca_pcs: int = 2,
    cca_compute_corrected_pvalues: bool = True,
    cca_coarse_start: float = 0.1,
    cca_coarse_end: float = 1.0,
    cca_coarse_step: float = 0.1,
    cca_fine_range: float = 0.02,
    cca_fine_step: float = 0.01,
) -> dict
```

### Parameters by stage

| Stage | Key parameters | Meaning |
| --- | --- | --- |
| Inputs | `atac_count_data_path`, `atac_sample_meta_path`, `atac_output_dir` | Input matrix/metadata and output root for ATAC branch |
| QC/filtering | `min_features`, `max_features`, `min_cells_per_sample` | ATAC feature/cell filtering |
| ATAC latent prep | `doublet_detection`, `tfidf_scale_factor`, `log_transform`, `drop_first_lsi`, `num_cell_hvfs` | ATAC-specific preprocessing |
| Cell types | `leiden_cluster_resolution`, `existing_cell_types`, `n_target_cell_clusters` | Cluster-based annotation |
| Sample embedding | `sample_hvg_number`, `sample_embedding_dimension`, `harmony_for_proportion` | Sample-level embedding settings |
| Resolution selection | `cca_*`, `trajectory_col`, `n_cca_pcs` | Optional CCA resolution search |

Returns `adata_cell`, `adata_sample`, `pseudo_adata`, `status_flags`.

## `preprocess(...)` and `preprocess_linux(...)`

- Source: `preparation/atac_preprocess_cpu.py` and `preparation/atac_preprocess_gpu.py`
- Purpose: preprocess scATAC data (feature QC, TF-IDF/LSI style workflow, optional doublet detection).

## `cell_types(...)` and `cell_types_linux(...)`

- Source: `preparation/cell_type_cpu.py` and `preparation/cell_type_gpu.py`
- Purpose: cell clustering and cell type labeling.

## `calculate_sample_embedding(..., atac=True)`

- Source: `sample_embedding/calculate_sample_embedding.py`
- Purpose: derive expression/proportion sample embeddings for ATAC branch.

## `find_optimal_cell_resolution(...)`

- Source: `parameter_selection/cpu_optimal_resolution.py` (GPU: `gpu_optimal_resolution.py`)
- Purpose: optimize cell type resolution by downstream CCA behavior.

## Related

- [RNA API](rna.md)
- [Multi-omics API](multiomics.md)
- [Downstream API](downstream.md)

