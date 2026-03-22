# Multi-omics API

## `multiomics_wrapper(...)`

```python
multiomics_wrapper(
    rna_file=None,
    atac_file=None,
    multiomics_output_dir=None,
    integration=True,
    integration_preprocessing=True,
    dimensionality_reduction=True,
    find_optimal_resolution=False,
    rna_sample_meta_file=None,
    atac_sample_meta_file=None,
    additional_hvg_file=None,
    rna_sample_column="sample",
    atac_sample_column="sample",
    sample_col="sample",
    batch_col=None,
    celltype_col="cell_type",
    modality_col="modality",
    multiomics_verbose=True,
    save_intermediate=True,
    use_gpu=True,
    random_state=42,
    run_glue_preprocessing=True,
    run_glue_training=True,
    run_glue_gene_activity=True,
    cell_type_cluster=True,
    run_glue_visualization=True,
    ensembl_release=98,
    species="homo_sapiens",
    use_highly_variable=True,
    n_top_genes=2000,
    n_pca_comps=50,
    n_lsi_comps=50,
    lsi_n_iter=15,
    gtf_by="gene_name",
    flavor="seurat_v3",
    generate_umap=False,
    compression="gzip",
    consistency_threshold=0.05,
    treat_sample_as_batch=True,
    save_prefix="glue",
    k_neighbors=10,
    use_rep="X_glue",
    metric="cosine",
    existing_cell_types=False,
    n_target_clusters=10,
    cluster_resolution=0.8,
    use_rep_celltype="X_glue",
    markers=None,
    generate_umap_celltype=True,
    plot_columns=None,
    min_cells_sample=1,
    min_cell_gene=10,
    min_features=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    doublet=True,
    sample_hvg_number=2000,
    preserve_cols_for_sample_embedding=None,
    n_expression_components=10,
    n_proportion_components=10,
    multiomics_harmony_for_proportion=True,
    optimization_target="rna",
    sev_col="sev.level",
    resolution_use_rep="X_glue",
    num_PCs=20,
    visualize_cell_types=True,
    coarse_start=0.1,
    coarse_end=1.0,
    coarse_step=0.1,
    fine_range=0.05,
    fine_step=0.01,
    analyze_modality_alignment=True,
    compute_corrected_pvalues=False,
    integrated_h5ad_path=None,
    pseudobulk_h5ad_path=None,
    status_flags=None,
) -> Dict[str, Any]
```

### Parameters by stage

| Stage | Key parameters | Meaning |
| --- | --- | --- |
| Inputs | `rna_file`, `atac_file`, `multiomics_output_dir` | Required RNA/ATAC inputs and output root |
| GLUE integration | `run_glue_preprocessing`, `run_glue_training`, `run_glue_gene_activity`, `run_glue_visualization` | Integration sub-steps |
| Integration preprocessing | `min_cells_sample`, `min_cell_gene`, `min_features`, `pct_mito_cutoff` | QC/filtering after integration |
| Cell type assignment | `cell_type_cluster`, `cluster_resolution`, `use_rep_celltype`, `generate_umap_celltype` | Multi-omics cell annotation |
| Sample embedding | `sample_hvg_number`, `n_expression_components`, `n_proportion_components` | Sample-level DR settings |
| Optimal resolution (API only) | `find_optimal_resolution`, `optimization_target`, `coarse_*`, `fine_*` | Resolution search for integrated data |

## `multiomics_preparation(...)`

- Source: `preparation/multi_omics_glue.py`
- Purpose: GLUE preprocessing/training/gene activity/visualization on RNA + ATAC inputs.

## `integrate_preprocess(...)`

- Source: `preparation/multi_omics_preprocess.py`
- Purpose: preprocess integrated object before sample-level analysis.

## `cell_types_multiomics(...)` and `cell_types_multiomics_linux(...)`

- Source: `preparation/multi_omics_cell_type_cpu.py` and `preparation/multi_omics_cell_type_gpu.py`
- Purpose: transfer/assign cell types in integrated multi-omics space.

## `calculate_multiomics_sample_embedding(...)`

- Source: `sample_embedding/calculate_multiomics_sample_embedding.py`
- Purpose: generate expression/proportion sample embeddings for integrated multi-omics data.

## `find_optimal_cell_resolution_multiomics(...)`

- Source: `parameter_selection/multi_omics_optimal_resolution_cpu.py` (GPU variant available)
- Purpose: multi-omics resolution search (documented here; not included in tutorials).

## `replace_optimal_dimension_reduction(...)`

- Source: `parameter_selection/multi_omics_unify_optimal.py`
- Purpose: swap optimized embedding outputs back into pseudobulk AnnData.

## Related

- [RNA API](rna.md)
- [ATAC API](atac.md)
- [Downstream API](downstream.md)

