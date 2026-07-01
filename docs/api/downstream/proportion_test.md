# `proportion_test`

Tests whether per-sample cell-type proportions differ across groups. Computes CLR (centered log-ratio) transformed proportions, then applies a limma-style empirical Bayes moderated test using `ebayes_test` (a Python port of the limma eBayes variance shrinkage). The CLR transform is compositional-aware — each sample is referenced against its own geometric mean over cell types, so passive shifts driven by the sum-to-one constraint do not inflate false positives. The test is run for every pairwise combination of groups, with BH-FDR applied globally across all (pair, cell-type) hypotheses. Groups can be supplied either via a column in `adata.obs` (`group_col`, preferred) or an explicit `sample_to_clade` mapping (from [`cluster`](cluster.md)).

**Source:** `sample_clustering/proportion_test.py:149`

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
| `celltype_col` | str | `"celltype"` | Cell-type column. Set to `"cell_type"` for SampleDisco pipelines. |
| `output_dir` | str, optional | `None` | Writes results + plots to `{output_dir}/proportion_test/`. |
| `verbose` | bool | `True` | Accepted for API symmetry; currently unused internally. |

## Returns

`dict` — `all_results`, a mapping `{"{group1}_vs_{group2}": DataFrame}`. Each DataFrame has columns `celltype`, `logFC`, `p_value`, `FDR` (sorted by `FDR`). When `output_dir` is given, the same results are also written to disk.

## Output files

Written directly under `{output_dir}` (created if missing):

- `proportion_test_{group1}_vs_{group2}.csv` — one file per pairwise comparison; columns `celltype`, `logFC`, `p_value`, `FDR`.
- `proportion_test_significant_summary.txt` — text summary of cell types with `FDR < 0.01`, per comparison.
- `proportion_heatmap_group_by_celltype.png` — heatmap of cell-type proportions per sample, grouped by group label.
- `proportion_heatmap_group_by_celltype_clr.png` — CLR-transformed version of the proportion heatmap.
- `proportion_heatmap_group_by_celltype_zscore.png` — z-scored version of the proportion heatmap.
- `proportion_significance_matrix.png` — group-pair × cell-type significance grid.
- `proportion_boxplot_{comparison}.png` — per-comparison cell-type proportion boxplots.

## Usage

```python
from sampledisco.sample_clustering.proportion_test import proportion_test

results = proportion_test(
    adata=adata_cell,
    sample_col="sample",
    sample_to_clade=clade_map,        # from cluster(...)
    celltype_col="cell_type",
    output_dir="sampledisco_demo_output/rna/sample_cluster",
)
```
