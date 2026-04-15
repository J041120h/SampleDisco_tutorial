# `raisinfit`

Python port of the RAISIN hierarchical generalized linear model for differential expression. `raisinfit` estimates mean expression and both cell-level and sample-level variance components given a sample × cell-type design, optionally correcting for batch with ComBat before fitting. The returned `fit` object is consumed by [`run_pairwise_tests`](run_pairwise_tests.md), which runs the actual pairwise contrasts and emits volcano plots. Supports unpaired, paired, continuous, and custom designs.

**Source:** `sample_clustering/RAISIN.py:244`

## Signature

```python
def raisinfit(
    adata,
    sample_col,
    testtype="unpaired",
    group_col=None,
    individual_col=None,
    batch_col=None,
    sample_to_clade=None,
    custom_design=None,
    intercept=True,
    filtergene=False,
    filtergenequantile=0.5,
    n_jobs=None,
    verbose=True,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Single-cell or pseudobulk AnnData. |
| `sample_col` | str | — | Sample identifier column. |
| `testtype` | str | `"unpaired"` | One of `"unpaired"`, `"paired"`, `"continuous"`, `"custom"`. |
| `group_col` | str, optional | `None` | Column with grouping/feature variable. Takes precedence over `sample_to_clade` when present. |
| `individual_col` | str, optional | `None` | Subject column for paired designs. |
| `batch_col` | str, optional | `None` | Triggers ComBat before fitting. |
| `sample_to_clade` | dict, optional | `None` | `{sample_id: group_label}` — used when `group_col` is absent. |
| `custom_design` | dict, optional | `None` | Required when `testtype="custom"`; keys `"X"`, `"Z"`, `"group"`. |
| `intercept` | bool | `True` | Include intercept in the fixed-effect design. |
| `filtergene` | bool | `False` | Drop lowly expressed genes before fitting. |
| `filtergenequantile` | float | `0.5` | Quantile threshold used when `filtergene=True`. |
| `n_jobs` | int, optional | `None` | CPU cores for parallel fitting (default: all cores). |
| `verbose` | bool | `True` | Print progress. |

## Returns

`dict` — the fit object with keys:

| Key | Meaning |
| --- | --- |
| `"mean"` | gene × sample expression mean matrix |
| `"sigma2"` | gene × group between-sample variance |
| `"omega2"` | gene × sample within-sample variance |
| `"X"` | fixed-effect design matrix |
| `"Z"` | random-effect design matrix |
| `"group"` | group assignment per sample |
| `"failgroup"` | groups where variance estimation failed |
| `"sample_names"` | ordered sample identifiers |
| `"batch_corrected"` | whether ComBat was applied |

## Usage

```python
from genodistance.sample_clustering import raisinfit, run_pairwise_tests

fit = raisinfit(
    adata=adata_cell,
    sample_col="sample",
    sample_to_clade=expr_clusters,
    testtype="unpaired",
    batch_col=None,
    intercept=True,
    n_jobs=8,
)

run_pairwise_tests(
    fit=fit,
    output_dir="/results/rna/raisin_results_expression",
    fdr_threshold=0.05,
)
```
