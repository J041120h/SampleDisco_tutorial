# RNA API

This page documents core RNA functions called by the pipeline.

## `preprocess(...)` (CPU)

Source: `code/preparation/rna_preprocess_cpu.py`

```python
preprocess(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
    cell_embedding_num_PCs=20,
    num_harmony_iterations=30,
    num_cell_hvgs=2000,
    min_cells=500,
    min_genes=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    cell_level_batch_key=None,
    verbose=True,
)
```

Preprocess RNA data and generate preprocessed `adata_cell` and `adata_sample`.

## `preprocess_linux(...)` (GPU)

Source: `code/preparation/rna_preprocess_gpu.py`

```python
preprocess_linux(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
    cell_embedding_num_PCs=20,
    num_harmony_iterations=30,
    num_cell_hvgs=2000,
    min_cells=500,
    min_genes=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    cell_level_batch_key=None,
    verbose=True,
)
```

GPU equivalent of RNA preprocessing.

## `cell_types(...)` (CPU)

Source: `code/preparation/cell_type_cpu.py`

```python
cell_types(
    anndata_cell,
    anndata_sample=None,
    cell_type_column="cell_type",
    existing_cell_types=False,
    n_target_clusters=None,
    umap=True,
    save=False,
    output_dir=None,
    defined_output_path=None,
    defined_sample_output_path=None,
    leiden_cluster_resolution=0.8,
    cell_embedding_column=None,
    cell_embedding_num_PCs=20,
    verbose=True,
    umap_plots=True,
    _recursion_depth=0,
)
```

Assign cell types from cell-level embedding and update sample-level object.

## `cell_types_linux(...)` (GPU)

Source: `code/preparation/cell_type_gpu.py`

```python
cell_types_linux(
    anndata_cell,
    anndata_sample=None,
    cell_type_column="cell_type",
    existing_cell_types=False,
    n_target_clusters=None,
    umap=True,
    save=False,
    output_dir=None,
    defined_output_path=None,
    defined_sample_output_path=None,
    leiden_cluster_resolution=0.8,
    cell_embedding_column=None,
    cell_embedding_num_PCs=20,
    verbose=True,
    umap_plots=True,
    _recursion_depth=0,
)
```

GPU equivalent of cell type assignment.

## `calculate_sample_embedding(...)`

Source: `code/sample_embedding/calculate_sample_embedding.py`

```python
calculate_sample_embedding(
    adata,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="./",
    sample_hvg_number=2000,
    combat_timeout=None,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    preserve_cols_in_sample_embedding=None,
    use_gpu=False,
    atac=False,
    save=True,
    verbose=True,
)
```

Create expression/proportion sample embeddings from preprocessed single-cell data.

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
    modality="rna",
    trajectory_col="sev.level",
    ...
)
```

Search for a clustering resolution maximizing CCA-based trajectory concordance.

