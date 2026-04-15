# `replace_optimal_dimension_reduction`

Folds the results of [`find_optimal_cell_resolution_multiomics_linux`](find_optimal_cell_resolution_multiomics_linux.md) back into the canonical pseudobulk object. Reads the per-DR-type optimal AnnData files and replaces the `X_DR_expression` / `X_DR_proportion` entries (plus the underlying pseudobulk expression matrix) in `pseudobulk_sample.h5ad`. After this step, downstream modules automatically pick up the optimal embeddings.

**Source:** `parameter_selection/multi_omics_unify_optimal.py:8`

## Signature

```python
def replace_optimal_dimension_reduction(
    base_path: str,
    expression_resolution_dir: Optional[str] = None,
    proportion_resolution_dir: Optional[str] = None,
    pseudobulk_path: Optional[str] = None,
    optimization_target: str = "rna",
    verbose: bool = True,
) -> sc.AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `base_path` | str | — | Multi-omics output root, e.g. `/results/multiomics`. |
| `expression_resolution_dir` | str, optional | `None` | Overrides `{base_path}/resolution_optimization_expression`. |
| `proportion_resolution_dir` | str, optional | `None` | Overrides `{base_path}/resolution_optimization_proportion`. |
| `pseudobulk_path` | str, optional | `None` | Overrides `{base_path}/pseudobulk/pseudobulk_sample.h5ad`. |
| `optimization_target` | str | `"rna"` | Target modality (`"rna"` or `"atac"`); used to construct expected file names such as `optimal_rna_expression.h5ad`. |
| `verbose` | bool | `True` | Print progress. |

## Returns

The updated pseudobulk AnnData (also written back to disk).

## Expected layout

```
{base_path}/
├── resolution_optimization_expression/
│   └── Integration_optimization_{target}_expression/summary/optimal_{target}_expression.h5ad
├── resolution_optimization_proportion/
│   └── Integration_optimization_{target}_proportion/summary/optimal_{target}_proportion.h5ad
└── pseudobulk/pseudobulk_sample.h5ad
```

## Usage

```python
from genodistance.parameter_selection import replace_optimal_dimension_reduction

updated_pseudo = replace_optimal_dimension_reduction(
    base_path="/results/multiomics",
    optimization_target="rna",
)
```
