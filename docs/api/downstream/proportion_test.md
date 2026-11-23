# `proportion_test`

Tests whether cell-type proportions differ between groups of samples. Per-sample proportions are computed from the cells, clipped away from 0 and 1, CLR-transformed (each sample against the geometric mean of its own proportions, which removes the sum-to-one coupling), and compared with a limma-style moderated t-test (`ebayes_test`, a Python port of limma's empirical-Bayes variance shrinkage). Every pair of groups is tested; BH-FDR is applied once across all (pair, cell type) tests. Pairs where either group has fewer than 2 samples are skipped with a warning.

Groups come from `group_col` (the most common value per sample) or, if `group_col` is `None`, from a `sample_to_clade` mapping such as the output of [`cluster`](cluster.md).

**Import:** `from sampledisco.sample_clustering.proportion_test import proportion_test`

**Source:** `sample_clustering/proportion_test.py:151`

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
| `adata` | AnnData | — | **Cell-level** AnnData with `sample_col` and `celltype_col` in `.obs`. |
| `sample_col` | str | — | Sample column. |
| `group_col` | str or None | `None` | Group column in `adata.obs`. When set, `sample_to_clade` is ignored (with a warning if both are given). Wrapper key: `*_cluster_differential_gene_group_col`. |
| `sample_to_clade` | dict or None | `None` | `{sample_id: group}`; samples not in it are left out. Used only when `group_col` is `None`. |
| `celltype_col` | str | `"celltype"` | Cell-type column. Note the default: SampleDisco's own labels are in `"cell_type"`. |
| `output_dir` | str or None | `None` | Write results and plots here (created if missing). |
| `verbose` | bool | `True` | Accepted; not used. |

At least 2 groups are required (`ValueError` otherwise).

## Returns

`dict`: `{"{group1}_vs_{group2}": DataFrame}` with columns `celltype`, `logFC` (CLR scale, group1 − group2), `p_value`, `FDR`, sorted by `FDR`.

## Output files

Directly under `output_dir`:

- `proportion_test_{group1}_vs_{group2}.csv`: one per tested pair.
- `proportion_test_significant_summary.txt`: cell types with **FDR < 0.01** per pair, and skipped pairs.
- `proportion_heatmap_group_by_celltype.png`: mean proportion per group × cell type.
- `proportion_heatmap_group_by_celltype_clr.png`, `proportion_heatmap_group_by_celltype_zscore.png`: the same on CLR and z-score scales, cell types ordered by significance.
- `proportion_significance_matrix.png`: pair × cell type, FDR < 0.01.
- `proportion_boxplot_{group1}_vs_{group2}.png`: per-sample proportions of the significant cell types, **only for pairs with at least one cell type at FDR < 0.01**.

## Usage

Continuing from the RNA chain, healthy (`sev.level` 1) vs severe (4):

```python
import anndata as ad
from sampledisco.sample_clustering.proportion_test import proportion_test

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")

res = proportion_test(
    adata,
    sample_col="sample",
    group_col="sev.level",
    celltype_col="cell_type",
    output_dir="api_demo/rna/proportion_test",
)
for pair, table in res.items():
    print(pair)
    print(table.round(3))
```

**Writes** → under `api_demo/rna/proportion_test/`: `proportion_test_1_vs_4.csv`, `proportion_test_significant_summary.txt`, `proportion_heatmap_group_by_celltype.png`, `proportion_heatmap_group_by_celltype_clr.png`, `proportion_heatmap_group_by_celltype_zscore.png`, `proportion_significance_matrix.png`. No cell type reaches FDR < 0.01, so no box plot is written. Printed (about 3 s):

```text
1_vs_4
  celltype  logFC  p_value    FDR
7        8 -2.335    0.005  0.041
6        7  1.839    0.011  0.043
4        5  1.490    0.022  0.059
1        2 -1.494    0.162  0.324
0        1 -0.435    0.460  0.614
5        6  0.473    0.416  0.614
3        4  0.292    0.552  0.631
2        3  0.170    0.745  0.745
```

`logFC` is group 1 (healthy) minus group 4 (severe) on the CLR scale.
