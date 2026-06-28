# `replace_optimal_dimension_reduction`

Folds the results of a per-DR-type resolution optimization back into the canonical pseudobulk object. Reads the optimal expression and proportion AnnData files, uses the optimal expression object as the template (carrying the updated `X`, `var`, and DR results), merges the original `obs` metadata plus non-DR `uns`/`obsm` from `pseudobulk_sample.h5ad`, overwrites the `X_DR_expression` / `X_DR_proportion` entries, writes a `.backup` of the original, and exports CSV embeddings under `{base_path}/embeddings/`.

!!! warning "Legacy helper — off the current pipeline path"
    This function still operates on the **legacy two-key** layout (`X_DR_expression` / `X_DR_proportion`). The current SampleDisco pipeline produces a **single** sample embedding, `adata.uns['X_DR_sample']`, written in place by [`compute_sample_embedding`](../shared/calculate_sample_embedding.md). It is not part of the config-driven wrapper path, and parameter selection is now handled by alpha/block-weight **autotune** (`sampledisco.parameter_selection.autotune.run_autotune`, enabled via `multiomics_autotune_enable` in the wrapper) rather than a resolution sweep. Keep this helper only if you are working with older two-key outputs.

**Source:** `utils/unify_optimal.py:7`

## Signature

```python
def replace_optimal_dimension_reduction(
    base_path: str,
    verbose: bool = True,
    modality: str = "RNA",
) -> sc.AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `base_path` | str | — | Run output directory, e.g. `.../covid_25_sample/rna`. Expects the optimization summary subdirs and `pseudobulk/pseudobulk_sample.h5ad` underneath (see layout below). |
| `verbose` | bool | `True` | Print progress. |
| `modality` | str | `"RNA"` | Target modality; must be `"RNA"` or `"ATAC"` (case-insensitive). Used to construct the `{MODALITY}_resolution_optimization_*` directory names. |

## Returns

The updated pseudobulk AnnData (also written back to `pseudobulk_sample.h5ad`, with the original preserved as a `.backup`).

## Expected layout

```
{base_path}/
├── {MODALITY}_resolution_optimization_expression/summary/optimal.h5ad
├── {MODALITY}_resolution_optimization_proportion/summary/optimal.h5ad
└── pseudobulk/pseudobulk_sample.h5ad
```

where `{MODALITY}` is `RNA` or `ATAC`. The three objects must share the same `obs` index order; the function raises `ValueError` if the sample indices do not match.

## Side effects

- Overwrites `pseudobulk_sample.h5ad` in place, after copying the original to `pseudobulk_sample.h5ad.backup` (skipped if a backup already exists).
- Writes `embeddings/sample_expression_embedding.csv` and `embeddings/sample_proportion_embedding.csv` under `base_path`.

## Usage

```python
from sampledisco.utils.unify_optimal import replace_optimal_dimension_reduction

updated_pseudo = replace_optimal_dimension_reduction(
    base_path="sampledisco_demo_output/multiomics/rna",
    modality="RNA",
)
```
