# `find_optimal_cell_resolution_multiomics_linux`

**Module**: `code/parameter_selection/multi_omics_optimal_resolution_gpu.py`

```python
find_optimal_cell_resolution_multiomics_linux(
    AnnData_integrated,
    output_dir,
    optimization_target="rna",
    dr_type="expression",
    sev_col="sev.level",
    batch_col=None,
    sample_col="sample",
    celltype_col="cell_type",
    modality_col="modality",
    use_rep="X_glue",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    preserve_cols=None,
    hvg_modality="RNA",
    coarse_start=0.1,
    coarse_end=1.0,
    coarse_step=0.1,
    fine_range=0.05,
    fine_step=0.01,
    visualize_cell_types=True,
    analyze_modality_alignment=True,
    compute_corrected_pvalues=False,
    verbose=True,
)
```

GPU/Linux two-pass resolution search for integrated multi-omics.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `AnnData_integrated` | required | Integrated RNA+ATAC AnnData. |
| `output_dir` | required | Output root for optimization artifacts. |
| `optimization_target` | `"rna"` | Modality score to optimize (`"rna"` or `"atac"`). |
| `dr_type` | `"expression"` | DR branch to optimize (`"expression"` or `"proportion"`). |
| `sev_col` | `"sev.level"` | Supervision column for CCA scoring. |
| `batch_col` | `None` | Batch covariate(s) for sample embedding recomputation. |
| `sample_col` | `"sample"` | Sample ID column. |
| `celltype_col` | `"cell_type"` | Cell type column. |
| `modality_col` | `"modality"` | Modality column. |
| `use_rep` | `"X_glue"` | Embedding used for clustering. |
| `sample_hvg_number` | `2000` | HVG count during sample embedding recomputation. |
| `n_expression_components` | `10` | Expression embedding dimensions. |
| `n_proportion_components` | `10` | Proportion embedding dimensions. |
| `harmony_for_proportion` | `True` | Harmony correction on proportion branch. |
| `preserve_cols` | `None` | Metadata columns to preserve. |
| `hvg_modality` | `"RNA"` | Modality used for HVG selection. |
| `coarse_start/end/step` | `0.1/1.0/0.1` | Coarse grid search range. |
| `fine_range/fine_step` | `0.05/0.01` | Fine search near best coarse resolution. |
| `visualize_cell_types` | `True` | Create visualizations in each resolution run. |
| `analyze_modality_alignment` | `True` | Compute RNA/ATAC alignment diagnostics. |
| `compute_corrected_pvalues` | `False` | Compute corrected p-values for selection bias. |
| `verbose` | `True` | Print detailed logs. |

## Returns

- `optimal_resolution`
- `results_df`

