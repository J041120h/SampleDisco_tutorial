# Multi-omics API

This page documents core multi-omics functions called by the pipeline.

## `multiomics_preparation(...)`

Source: `code/preparation/multi_omics_glue.py`

```python
multiomics_preparation(
    rna_file,
    atac_file,
    rna_sample_meta_file=None,
    atac_sample_meta_file=None,
    additional_hvg_file=None,
    run_preprocessing=True,
    run_training=True,
    run_gene_activity=True,
    run_visualization=True,
    ...
)
```

Run GLUE preprocessing, training, and optional gene activity / visualization.

## `integrate_preprocess(...)`

Source: `code/preparation/multi_omics_preprocess.py`

```python
integrate_preprocess(
    output_dir,
    h5ad_path=None,
    sample_column="sample",
    modality_col="modality",
    min_cells_sample=1,
    min_cell_gene=10,
    min_features=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    doublet=True,
    ...
)
```

Preprocess integrated RNA+ATAC object for sample-level embedding.

## `cell_types_multiomics(...)` (CPU)

Source: `code/preparation/multi_omics_cell_type_cpu.py`

```python
cell_types_multiomics(
    adata,
    modality_column="modality",
    rna_modality_value="RNA",
    atac_modality_value="ATAC",
    cell_type_column="cell_type",
    cluster_resolution=0.8,
    use_rep="X_glue",
    ...
)
```

Assign cell types on integrated embeddings.

## `cell_types_multiomics_linux(...)` (GPU)

Source: `code/preparation/multi_omics_cell_type_gpu.py`

```python
cell_types_multiomics_linux(
    adata,
    modality_column="modality",
    rna_modality_value="RNA",
    atac_modality_value="ATAC",
    cell_type_column="cell_type",
    cluster_resolution=0.8,
    use_rep="X_glue",
    ...
)
```

GPU equivalent of multi-omics cell type assignment.

## `calculate_multiomics_sample_embedding(...)`

Source: `code/sample_embedding/calculate_multiomics_sample_embedding.py`

```python
calculate_multiomics_sample_embedding(
    adata,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="./",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    hvg_modality="RNA",
    modality_col="modality",
    ...
)
```

Create expression/proportion sample embeddings for integrated multi-omics data.

## `find_optimal_cell_resolution_multiomics(...)` and `find_optimal_cell_resolution_multiomics_linux(...)`

Sources:

- `code/parameter_selection/multi_omics_optimal_resolution_cpu.py`
- `code/parameter_selection/multi_omics_optimal_resolution_gpu.py`

```python
find_optimal_cell_resolution_multiomics(
    AnnData_integrated,
    output_dir,
    optimization_target="rna",
    dr_type="expression",
    sev_col="sev.level",
    ...
)
```

Search optimal multi-omics cell type resolution by CCA-driven criteria.

## `replace_optimal_dimension_reduction(...)`

Source: `code/parameter_selection/multi_omics_unify_optimal.py`

```python
replace_optimal_dimension_reduction(
    base_path,
    expression_resolution_dir=None,
    proportion_resolution_dir=None,
    pseudobulk_path=None,
    optimization_target="rna",
    verbose=True,
)
```

Replace default pseudobulk DR results with selected optimal-resolution outputs.

