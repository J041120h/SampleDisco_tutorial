# `run_trajectory_gam_differential_gene_analysis`

Fits a Generalized Additive Model (GAM) per gene with pseudotime as the smooth predictor and optional covariates as fixed effects. After fitting, effect sizes (max minus min of the smoothed curve) and F-test p-values are computed per gene, corrected with BH FDR, and filtered by `fdr_threshold` and `effect_size_threshold`. Lamian-style visualizations are generated: heatmap, volcano, per-gene curves, per-sample curves, cluster patterns, and a one-page results summary.

**Source:** `sample_trajectory/trajectory_diff_gene.py:758`

## Signature

```python
def run_trajectory_gam_differential_gene_analysis(
    pseudobulk_adata: ad.AnnData,
    pseudotime_source: Union[str, pd.DataFrame, Dict],
    *,
    sample_col: str = "sample",
    pseudotime_col: str = "pseudotime",
    covariate_columns: Optional[List[str]] = None,
    fdr_threshold: float = 0.01,
    effect_size_threshold: float = 1.0,
    top_n_genes: int = 100,
    num_splines: int = 5,
    spline_order: int = 3,
    output_dir: str = "trajectory_diff_gene_results_single",
    visualization_gene_list: Optional[List[str]] = None,
    generate_visualizations: bool = True,
    group_col: Optional[str] = None,
    n_clusters: int = 3,
    top_n_genes_for_curves: int = 20,
    verbose: bool = True,
) -> pd.DataFrame
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudobulk_adata` | AnnData | — | Sample × gene pseudobulk. |
| `pseudotime_source` | str / DataFrame / dict | — | Path to CSV/TSV, DataFrame, or `{sample_id: pseudotime}` dict. |
| `sample_col` | str | `"sample"` | Sample identifier column in `pseudobulk_adata.obs`. |
| `pseudotime_col` | str | `"pseudotime"` | Column within `pseudotime_source` holding the values. |
| `covariate_columns` | list, optional | `None` | Additional fixed-effect covariates (batch, age, sex...). |
| `fdr_threshold` | float | `0.01` | FDR cutoff for pseudoDEG selection. |
| `effect_size_threshold` | float | `1.0` | Minimum max-smooth difference along pseudotime. |
| `top_n_genes` | int | `100` | Retain at most this many top-ranked pseudoDEGs. |
| `num_splines` | int | `5` | GAM basis size. |
| `spline_order` | int | `3` | Spline polynomial order (`2` for RNA in the canonical config, `3` for ATAC). |
| `output_dir` | str | — | Base directory for the results. |
| `visualization_gene_list` | list, optional | `None` | Named genes to always visualize, regardless of rank. |
| `generate_visualizations` | bool | `True` | Turn off to skip the visualization bundle. |
| `group_col` | str, optional | `None` | Optional metadata column for cluster comparisons in visualizations. |
| `n_clusters` | int | `3` | Number of gene clusters in the heatmap. |
| `top_n_genes_for_curves` | int | `20` | Genes shown in per-gene curve plots. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`pd.DataFrame` — per-gene result table with columns including `gene`, `effect_size`, `pvalue`, `FDR`, `rank`, `is_pseudoDEG`.

## Output files

Under `{output_dir}/`:

- `pseudoDEGs.csv` — filtered pseudoDEG table.
- `all_gene_results.csv` — full ranking.
- `visualizations/tde_heatmap.png`
- `visualizations/gene_curves.png`
- `visualizations/volcano_plot.png`
- `visualizations/sample_density.png`
- `visualizations/sample_level_curves.png`
- `visualizations/cluster_patterns.png`
- `visualizations/results_summary.png`

## Usage

```python
from genodistance.sample_trajectory import run_trajectory_gam_differential_gene_analysis

results_df = run_trajectory_gam_differential_gene_analysis(
    pseudobulk_adata=pseudo_adata,
    pseudotime_source=expr_pseudotime_dict,
    sample_col="sample",
    pseudotime_col="pseudotime",
    fdr_threshold=0.05,
    effect_size_threshold=1.0,
    top_n_genes=100,
    num_splines=5,
    spline_order=2,
    output_dir="/results/rna/trajectoryDEG/expression",
)
```
