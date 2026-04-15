# `proportion_test`

Tests whether per-sample cell-type proportions differ across groups. Computes logit-transformed proportions, then applies a limma-style empirical Bayes moderated F-test using `ebayes_test` (a Python port of the limma eBayes variance shrinkage). Groups can be supplied either via a column in `adata.obs` (`group_col`, preferred) or an explicit `sample_to_clade` mapping (from [`cluster`](cluster.md)).

**Source:** `sample_clustering/proportion_test.py:175`

## Signature

```python
def proportion_test(
    adata,
    sample_col,
    group_col=None,
    sample_to_clade=None,
    celltype_col="celltype",
    output_dir=None,
    verbose=True,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level AnnData whose `.obs` contains `sample_col` and `celltype_col`. |
| `sample_col` | str | — | Sample identifier column. |
| `group_col` | str, optional | `None` | Column in `adata.obs` defining the group label. Takes precedence over `sample_to_clade`. |
| `sample_to_clade` | dict, optional | `None` | `{sample_id: group_label}` mapping. Used only when `group_col` is `None`. |
| `celltype_col` | str | `"celltype"` | Cell-type column. Set to `"cell_type"` for SampleDisc pipelines. |
| `output_dir` | str, optional | `None` | Writes results + plots to `{output_dir}/proportion_test/`. |
| `verbose` | bool | `True` | Accepted for API symmetry; currently unused internally. |

## Returns

`None`. Side effect: writes results to disk.

## Output files

Under `{output_dir}/proportion_test/`:

- `proportion_test_results.csv` — per cell-type test statistics (F, p, FDR, effect size).
- `proportion_heatmap_group_by_celltype.png` — heatmap of cell-type proportions per sample, grouped by cluster.
- `proportion_significance_matrix.png` — cluster × cell-type significance grid.

## Usage

```python
from genodistance.sample_clustering import proportion_test

proportion_test(
    adata=adata_cell,
    sample_col="sample",
    sample_to_clade=expr_clusters,        # from cluster(...)
    celltype_col="cell_type",
    output_dir="/results/rna/sample_cluster/expression",
)
```
