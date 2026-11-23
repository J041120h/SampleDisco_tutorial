# `replace_optimal_dimension_reduction`

!!! danger "Legacy: not part of the current pipeline"
    This helper belongs to the pre-0.3 **two-key** layout (`X_DR_expression` / `X_DR_proportion`, a `pseudobulk/pseudobulk_sample.h5ad` file and a cell-resolution sweep). Nothing in SampleDisco 0.5 calls it, and none of the current outputs have the files it expects. The current pipeline writes one sample embedding, `uns['X_DR_sample']`, with [`compute_sample_embedding`](../shared/calculate_sample_embedding.md), and tunes the RMD weight with [`run_autotune`](../shared/run_autotune.md) instead of sweeping resolutions. Use this function only to post-process output folders written by those old releases. It may be removed in a future release.

Folds the results of an old per-DR-type resolution optimization back into the pseudobulk object. It reads the optimal expression and proportion AnnData files, uses the optimal expression object as the template (its `X`, `var` and DR results), merges in the original `obs` and the non-DR `uns`/`obsm` entries of `pseudobulk_sample.h5ad`, overwrites `X_DR_expression` / `X_DR_proportion`, keeps a `.backup` of the original, and exports CSV embeddings under `{base_path}/embeddings/`.

**Import:** `from sampledisco.utils.unify_optimal import replace_optimal_dimension_reduction`

**Source:** `utils/unify_optimal.py:7`

## Signature

```python
def replace_optimal_dimension_reduction(
    base_path: str,
    verbose: bool = True,
    modality: str = "RNA",
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `base_path` | str | — | An old run's output directory for one modality, e.g. `.../covid_25_sample/rna`, laid out as below. |
| `verbose` | bool | `True` | Print progress. |
| `modality` | str | `"RNA"` | `"RNA"` or `"ATAC"` (case-insensitive); used to build the `{MODALITY}_resolution_optimization_*` folder names. |

## Expected layout

```
{base_path}/
├── {MODALITY}_resolution_optimization_expression/summary/optimal.h5ad
├── {MODALITY}_resolution_optimization_proportion/summary/optimal.h5ad
└── pseudobulk/pseudobulk_sample.h5ad
```

The three objects must share the same `obs` index; otherwise a `ValueError` is raised.

## Returns

The updated pseudobulk `AnnData`, also written back to `pseudobulk_sample.h5ad`.

## Side effects

- Overwrites `pseudobulk/pseudobulk_sample.h5ad`, after copying the original to `pseudobulk_sample.h5ad.backup` (skipped if a backup already exists).
- Writes `embeddings/sample_expression_embedding.csv` and `embeddings/sample_proportion_embedding.csv` under `base_path`.

## Usage

There is no runnable example: the demo data and every 0.5 output lack the layout above. On an old output folder the call is `replace_optimal_dimension_reduction(base_path="<old run>/rna", modality="RNA")`.
