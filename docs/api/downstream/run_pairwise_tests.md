# `run_pairwise_tests`

Runs all pairwise RAISIN contrasts for the groups in a fitted model and generates the volcano plots, per-pair results CSVs, and summary figures. Internally calls `raisintest` once per pair (or once per non-control group when `control_group` is set) and aggregates the outputs into a `summary_plots/` folder with a cross-cluster pseudobulk heatmap and a summary dotplot of DE gene counts.

**Source:** `sampledisco/sample_clustering/RAISIN_TEST.py:609`

## Signature

```python
def run_pairwise_tests(
    fit,
    output_dir,
    groups_to_compare=None,
    control_group=None,
    fdrmethod="fdr_bh",
    n_permutations=100,
    fdr_threshold=0.05,
    top_n_genes: int = 50,
    make_summary_plots: bool = True,
    verbose=True,
    random_state: int = 42,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `fit` | dict | — | Output of [`raisinfit`](raisinfit.md). |
| `output_dir` | str | — | Writes to `{output_dir}/{group1}_vs_{group2}/` for each pair plus `{output_dir}/summary_plots/`. |
| `groups_to_compare` | list, optional | `None` | Subset of groups to consider. Defaults to all groups in `fit["group"]`. |
| `control_group` | str, optional | `None` | If set, tests every other group against this reference rather than all pairwise combinations. |
| `fdrmethod` | str | `"fdr_bh"` | statsmodels multi-testing method. |
| `n_permutations` | int | `100` | Number of permutations used by `raisintest` for empirical calibration. |
| `fdr_threshold` | float | `0.05` | Cutoff used when drawing volcano significance bands and counting DE genes. |
| `top_n_genes` | int | `50` | How many top genes to annotate in the volcano and heatmap. |
| `make_summary_plots` | bool | `True` | Emit `summary_plots/pseudobulk_heatmap.png`, `summary_plots/summary_dotplot.png`, `summary_plots/cluster_gene_zscore.png`. |
| `verbose` | bool | `True` | Print progress. |
| `random_state` | int | `42` | Seed passed to `raisintest` for the permutation calibration of each pair. |

## Returns

`(results_summary, all_results)` — a 2-tuple. `results_summary` is a `dict` mapping each `{group1}_vs_{group2}` comparison to its count of significant genes (`FDR < fdr_threshold`); `all_results` is a `dict` mapping each comparison to its per-pair RAISIN results `DataFrame`. The function also writes the per-pair CSVs/PNGs and the summary folder described below.

## Output files

Per pair:

- `{group1}_vs_{group2}/raisin_results.csv`
- `{group1}_vs_{group2}/volcano.png`
- `{group1}_vs_{group2}/volcano_labeled.png`

Across all pairs (a `summary_plots/` folder is always written; the PNG plots require `make_summary_plots=True`):

- `summary_plots/raisin_summary.csv`
- `summary_plots/raisin_summary.txt`
- `summary_plots/pseudobulk_heatmap.png`
- `summary_plots/summary_dotplot.png`
- `summary_plots/cluster_gene_zscore.png`
- `summary_plots/all_results_combined.csv`

## Usage

```python
from sampledisco.sample_clustering.RAISIN_TEST import run_pairwise_tests

results_summary, all_results = run_pairwise_tests(
    fit=fit,
    output_dir="sampledisco_demo_output/rna/raisin_results",
    groups_to_compare=None,
    control_group=None,
    fdr_threshold=0.05,
    top_n_genes=50,
    make_summary_plots=True,
)
```
