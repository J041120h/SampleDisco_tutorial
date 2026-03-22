# ATAC API

This page documents core ATAC functions called by the pipeline.

## `preprocess(...)` (CPU)

Source: `code/preparation/atac_preprocess_cpu.py`

```python
preprocess(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
    cell_embedding_num_PCs=50,
    num_harmony_iterations=30,
    num_cell_hvfs=50000,
    min_cells=1,
    min_features=2000,
    max_features=15000,
    min_cells_per_sample=1,
    exclude_features=None,
    cell_level_batch_key=None,
    doublet_detection=True,
    tfidf_scale_factor=1e4,
    log_transform=True,
    drop_first_lsi=True,
    verbose=True,
)
```

Preprocess ATAC data, including feature QC and LSI-ready representation.

## `preprocess_linux(...)` (GPU)

Source: `code/preparation/atac_preprocess_gpu.py`

```python
preprocess_linux(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
    cell_embedding_num_PCs=50,
    num_harmony_iterations=30,
    num_cell_hvfs=50000,
    min_cells=1,
    min_features=2000,
    max_features=15000,
    min_cells_per_sample=1,
    exclude_features=None,
    cell_level_batch_key=None,
    doublet_detection=True,
    tfidf_scale_factor=1e4,
    log_transform=True,
    drop_first_lsi=True,
    verbose=True,
)
```

GPU equivalent of ATAC preprocessing.

## `cell_types(...)` and `cell_types_linux(...)`

Sources:

- `code/preparation/cell_type_cpu.py`
- `code/preparation/cell_type_gpu.py`

Same function family used in RNA, applied to ATAC preprocessed embeddings.

## `calculate_sample_embedding(...)` with `atac=True`

Source: `code/sample_embedding/calculate_sample_embedding.py`

```python
calculate_sample_embedding(
    adata,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="./",
    sample_hvg_number=50000,
    n_expression_components=30,
    n_proportion_components=30,
    harmony_for_proportion=True,
    preserve_cols_in_sample_embedding=None,
    use_gpu=False,
    atac=True,
    save=True,
    verbose=True,
)
```

Build sample-level expression/proportion embeddings for ATAC.

## `find_optimal_cell_resolution(...)` and `find_optimal_cell_resolution_linux(...)`

Sources:

- `code/parameter_selection/cpu_optimal_resolution.py`
- `code/parameter_selection/gpu_optimal_resolution.py`

```python
find_optimal_cell_resolution(
    adata_cell,
    adata_sample,
    output_dir,
    column,
    modality="atac",
    trajectory_col="sev.level",
    ...
)
```

Optimize cell type resolution using CCA-based downstream signal.

