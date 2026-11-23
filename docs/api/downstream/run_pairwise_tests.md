# `run_pairwise_tests`

Tests every pair of groups in a [`raisinfit`](raisinfit.md) fit (or every group against `control_group`) with `raisintest`: a contrast of the two group coefficients, a test statistic per gene, and p-values calibrated against `n_permutations` permuted statistics (t or normal, whichever fits the permutations better), then FDR. Each comparison gets its own folder with a results table and volcano plots; summary tables and plots across comparisons follow. A comparison whose group is not in the design, or whose test fails, is skipped and marked with a `FAILED` file.

!!! note "Groups with a single sample (0.6.0)"
    `raisinfit` cannot estimate the variance of a group with no residual degrees of freedom (e.g. a single sample). It lists such groups in `fit['failgroup']` and sets their variance to 0, which would make every gene look significant. In 0.6.0, `run_pairwise_tests` warns once (`RuntimeWarning`), skips every comparison whose variance uses such a group, writes `{g1}_vs_{g2}/FAILED` ("variance not estimable … not tested"), and records it in `raisin_summary.csv` with an empty `n_sig`. Skipped comparisons are not in the returned dicts. Calling `raisintest` directly on such a contrast returns `stat`, `pvalue` and `FDR` as `NaN` (the `Foldchange` is kept) and writes `raisin_results.csv` but no volcano plot. With k-means groups from [`cluster`](cluster.md), keep k low enough that every cluster has at least 2 samples.

**Import:** `from sampledisco.sample_clustering.RAISIN_TEST import run_pairwise_tests`

**Source:** `sample_clustering/RAISIN_TEST.py:642`

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
| `output_dir` | str | — | Root folder (see Output files). |
| `groups_to_compare` | list or None | `None` | Groups to include; `None` = all groups of the fit. |
| `control_group` | str or None | `None` | Test each other group against this one instead of all pairs. |
| `fdrmethod` | str | `"fdr_bh"` | `statsmodels` multiple-testing method. |
| `n_permutations` | int | `100` | Permutations for the null calibration of each comparison. |
| `fdr_threshold` | float | `0.05` | Significance cutoff for counts, volcano colors and summaries. |
| `top_n_genes` | int | `50` | Genes listed per comparison in `raisin_summary.txt` and shown in the heatmap and z-score plot. The labelled volcano always labels 5 genes. |
| `make_summary_plots` | bool | `True` | Write the summary plots and `all_results_combined.csv`. |
| `verbose` | bool | `True` | Print progress. |
| `random_state` | int | `42` | Seed of the permutations. |

## Returns

`(results_summary, all_results)`:

- `results_summary`: `{"{g1}_vs_{g2}": number of genes with FDR < fdr_threshold}`.
- `all_results`: `{"{g1}_vs_{g2}": DataFrame}` indexed by gene, columns `Foldchange`, `stat`, `pvalue`, `FDR`, sorted by FDR.

## Output files

Per comparison, under `{output_dir}/{g1}_vs_{g2}/`: `raisin_results.csv`, `volcano.png`, `volcano_labeled.png` (or `FAILED`).

Under `{output_dir}/summary_plots/`:

- Always (when at least one comparison was attempted): `raisin_summary.csv` (`comparison`, `status`, `n_genes`, `n_sig`, `n_sig_up`, `n_sig_down`) and `raisin_summary.txt`.
- With `make_summary_plots=True`: `cluster_gene_zscore.png` (falls back to the top-ranked genes when none is significant), `all_results_combined.csv`, and, only when some gene passes `fdr_threshold`, `pseudobulk_heatmap.png` and `summary_dotplot.png`.

## Usage

Continuing from the RNA chain (the fit is repeated so the example runs on its own):

```python
import anndata as ad
from sampledisco.sample_clustering.RAISIN import raisinfit
from sampledisco.sample_clustering.RAISIN_TEST import run_pairwise_tests

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
fit = raisinfit(adata, sample_col="sample", group_col="sev.level", n_jobs=8, verbose=False)

summary, results = run_pairwise_tests(fit, "api_demo/rna/raisin_results")
print(summary)
for pair, table in results.items():
    print(pair)
    print(table.head())
```

**Writes** → `api_demo/rna/raisin_results/1_vs_4/raisin_results.csv`, `volcano.png`, `volcano_labeled.png`, and under `summary_plots/`: `raisin_summary.csv`, `raisin_summary.txt`, `cluster_gene_zscore.png`, `all_results_combined.csv`. No gene reaches FDR < 0.05, so the heatmap and dot plot are skipped with a warning. Printed:

```text
{'1_vs_4': np.int64(0)}
1_vs_4
        Foldchange      stat    pvalue       FDR
IGKC     -1.394216 -4.174512  0.000030  0.296715
ZZEF1    -0.023033 -0.149775  0.880942  0.999956
ZYX      -0.137261 -0.727488  0.466927  0.999956
ZYG11B   -0.038858 -0.405257  0.685288  0.999956
ZXDC      0.024459  0.254218  0.799327  0.999956
```

About 80 s (both calls) and 4.6 GB including workers. Four samples per group leave RAISIN little power; the [RAISIN tutorial](../../tutorials/downstream/raisin_dge.md) and the [Sound Life walkthrough](../../tutorials/soundlife_rna.md) cover larger designs.
