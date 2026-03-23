# `run_pairwise_tests`

**Module**: `code/sample_clustering/RAISIN_TEST.py`

```python
run_pairwise_tests(
    fit,
    output_dir,
    groups_to_compare=None,
    control_group=None,
    fdrmethod="fdr_bh",
    n_permutations=100,
    fdr_threshold=0.05,
    top_n_genes=50,
    make_summary_plots=True,
    verbose=True,
)
```

Runs pairwise RAISIN comparisons and generates summary outputs.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `fit` | required | Output object from `raisinfit`. |
| `output_dir` | required | Directory for pairwise results and plots. |
| `groups_to_compare` | `None` | Group list to include; defaults to all groups in fit. |
| `control_group` | `None` | If set, compares each group against this control only. |
| `fdrmethod` | `"fdr_bh"` | Multiple-testing correction method. |
| `n_permutations` | `100` | Permutation count in RAISIN test step. |
| `fdr_threshold` | `0.05` | Significance threshold for summaries. |
| `top_n_genes` | `50` | Top genes per comparison for summary visualizations. |
| `make_summary_plots` | `True` | Generate heatmap/dotplot summary assets. |
| `verbose` | `True` | Print progress logs. |

## Returns

- `results_summary`
- `all_results`

