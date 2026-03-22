# RNA API

## `rna_wrapper(...)`

```python
rna_wrapper(
    rna_count_data_path: str = None,
    rna_output_dir: str = None,
    preprocessing: bool = True,
    cell_type_cluster: bool = True,
    derive_sample_embedding: bool = True,
    cca_based_cell_resolution_selection: bool = False,
    use_gpu: bool = False,
    verbose: bool = True,
    status_flags: dict = None,
    adata_cell_path: str = None,
    adata_sample_path: str = None,
    rna_sample_meta_path: str = None,
    cell_meta_path: str = None,
    pseudo_adata_path: str = None,
    sample_col: str = "sample",
    sample_level_batch_col: Optional[List[str]] = None,
    celltype_col: str = "cell_type",
    min_cells: int = 500,
    min_genes: int = 500,
    pct_mito_cutoff: float = 20,
    exclude_genes: list = None,
    num_cell_hvgs: int = 2000,
    cell_embedding_num_pcs: int = 20,
    num_harmony_iterations: int = 30,
    cell_level_batch_key: list = None,
    leiden_cluster_resolution: float = 0.8,
    cell_embedding_column: str = None,
    existing_cell_types: bool = False,
    n_target_cell_clusters: int = None,
    umap: bool = False,
    sample_hvg_number: int = 2000,
    sample_embedding_dimension: int = 10,
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
| Inputs | `rna_count_data_path`, `rna_sample_meta_path`, `rna_output_dir` | Input matrix/metadata and output root for RNA branch |
| Preprocessing | `min_cells`, `min_genes`, `pct_mito_cutoff`, `num_cell_hvgs` | QC and feature-selection settings |
| Batch handling | `sample_level_batch_col`, `cell_level_batch_key`, `num_harmony_iterations` | Batch columns and correction iterations |
| Cell types | `leiden_cluster_resolution`, `existing_cell_types`, `n_target_cell_clusters` | Cell label inference/reuse |
| Sample embedding | `sample_hvg_number`, `sample_embedding_dimension`, `harmony_for_proportion` | Sample-level DR settings |
| Resolution selection | `cca_*`, `n_cca_pcs`, `trajectory_col` | Optional CCA-based optimal resolution search |

Returns `adata_cell`, `adata_sample`, `pseudo_adata`, `status_flags`.

## `preprocess(...)` and `preprocess_linux(...)`

- Source: `preparation/rna_preprocess_cpu.py` and `preparation/rna_preprocess_gpu.py`
- Purpose: preprocess cell-level RNA data and generate preprocessed `adata_cell`/`adata_sample`.
- Typical config keys: `rna_min_cells`, `rna_min_genes`, `rna_pct_mito_cutoff`, `rna_num_cell_hvgs`.

## `cell_types(...)` and `cell_types_linux(...)`

- Source: `preparation/cell_type_cpu.py` and `preparation/cell_type_gpu.py`
- Purpose: cell clustering/cell type assignment on preprocessed data.
- Typical config keys: `rna_leiden_cluster_resolution`, `rna_existing_cell_types`, `rna_umap`.

## `calculate_sample_embedding(...)`

- Source: `sample_embedding/calculate_sample_embedding.py`
- Purpose: create pseudobulk expression/proportion representations and sample embeddings.
- Typical config keys: `rna_sample_hvg_number`, `rna_sample_embedding_dimension`, `rna_harmony_for_proportion`.

## `find_optimal_cell_resolution(...)`

- Source: `parameter_selection/cpu_optimal_resolution.py` (GPU: `gpu_optimal_resolution.py`)
- Purpose: resolution search by maximizing trajectory-related CCA signal.
- Typical config keys: `rna_cca_coarse_*`, `rna_cca_fine_*`, `rna_cca_compute_corrected_pvalues`.

## Related

- [ATAC API](atac.md)
- [Multi-omics API](multiomics.md)
- [Downstream API](downstream.md)

