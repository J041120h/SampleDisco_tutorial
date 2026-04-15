# `find_optimal_cell_resolution_multiomics_linux`

Two-pass resolution search for the integrated RNA + ATAC object. For a chosen `optimization_target` (RNA or ATAC) and `dr_type` (expression or proportion), the function sweeps Leiden resolutions on the joint embedding, rebuilds the multi-omics sample embedding at each resolution, and picks the value that maximizes CCA correlation with the trajectory column. GPU-accelerated throughout via [`cell_types_multiomics_linux`](cell_types_multiomics_linux.md) and [`calculate_multiomics_sample_embedding`](calculate_multiomics_sample_embedding.md).

**Source:** `parameter_selection/multi_omics_optimal_resolution_gpu.py:543`

## Signature

```python
def find_optimal_cell_resolution_multiomics_linux(
    AnnData_integrated: AnnData,
    output_dir: str,
    optimization_target: Literal["rna", "atac"] = "rna",
    dr_type: Literal["expression", "proportion"] = "expression",
    sev_col: str = "sev.level",
    batch_col: Optional[Union[str, List[str]]] = None,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    modality_col: str = "modality",
    use_rep: str = "X_glue",
    sample_hvg_number: int = 2000,
    n_expression_components: int = 10,
    n_proportion_components: int = 10,
    harmony_for_proportion: bool = True,
    preserve_cols: Optional[Union[str, List[str]]] = None,
    hvg_modality: str = "RNA",
    coarse_start: float = 0.1,
    coarse_end: float = 1.0,
    coarse_step: float = 0.1,
    fine_range: float = 0.05,
    fine_step: float = 0.01,
    visualize_cell_types: bool = True,
    analyze_modality_alignment: bool = True,
    compute_corrected_pvalues: bool = False,
    verbose: bool = True,
) -> Tuple[float, pd.DataFrame]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `AnnData_integrated` | AnnData | — | Integrated RNA + ATAC object with the joint embedding in `use_rep`. |
| `output_dir` | str | — | Results land under `{output_dir}/resolution_optimization_{dr_type}/`. |
| `optimization_target` | `"rna"` or `"atac"` | `"rna"` | Which modality's phenotype is used for CCA. |
| `dr_type` | `"expression"` or `"proportion"` | `"expression"` | Which sample embedding to optimize. |
| `sev_col` | str | `"sev.level"` | Trajectory/phenotype column. |
| `batch_col` | str or list, optional | `None` | Batch column(s) for sample embedding. |
| `sample_col` | str | `"sample"` | Sample identifier. |
| `celltype_col` | str | `"cell_type"` | Cell-type label column (rewritten per resolution). |
| `modality_col` | str | `"modality"` | Column identifying modality. |
| `use_rep` | str | `"X_glue"` | Representation for clustering. |
| `sample_hvg_number` | int | `2000` | HVGs per cell type in pseudobulk. |
| `n_expression_components` / `n_proportion_components` | int | `10 / 10` | DR components. |
| `harmony_for_proportion` | bool | `True` | Harmony on the proportion embedding. |
| `preserve_cols` | str or list, optional | `None` | Metadata carried into pseudobulk. |
| `hvg_modality` | str | `"RNA"` | Which modality's HVGs feed the pseudobulk. |
| `coarse_start` / `coarse_end` / `coarse_step` | float | `0.1 / 1.0 / 0.1` | Coarse grid. |
| `fine_range` / `fine_step` | float | `0.05 / 0.01` | Fine grid around the coarse winner. |
| `visualize_cell_types` | bool | `True` | Save per-resolution cell-type UMAPs. |
| `analyze_modality_alignment` | bool | `True` | Report modality alignment quality. |
| `compute_corrected_pvalues` | bool | `False` | Permutation-corrected p-values. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`Tuple[float, pd.DataFrame]` — `(best_resolution, trials_df)` with per-resolution CCA scores, p-values, and diagnostics.

## Output files

- `{output_dir}/resolution_optimization_{expression|proportion}/Integration_optimization_{target}_{dr_type}/` — per-resolution artifacts.
- `.../summary/optimal_{target}_{dr_type}.h5ad` — AnnData at the winning resolution.
- `.../summary/trials.csv`.

## Usage

```python
from genodistance.parameter_selection import find_optimal_cell_resolution_multiomics_linux

best_res, trials = find_optimal_cell_resolution_multiomics_linux(
    AnnData_integrated=adata_integrated,
    output_dir="/results/multiomics",
    optimization_target="rna",
    dr_type="expression",
    sev_col="sev.level",
    use_rep="X_glue",
)
```

After running, call [`replace_optimal_dimension_reduction`](replace_optimal_dimension_reduction.md) to fold the winning embedding back into `pseudobulk_sample.h5ad`.
