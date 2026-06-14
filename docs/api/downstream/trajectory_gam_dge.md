# `run_trajectory_gam_differential_gene_analysis`

Takes a **cell-level** `AnnData`, builds a per-sample (optionally per-cell-type) pseudobulk internally, then fits a Generalized Additive Model (GAM) per gene with pseudotime as the smooth predictor and optional covariates as fixed effects. After fitting, effect sizes and deviance/significance statistics are computed per gene, corrected with BH FDR, and used to select pseudoDEGs by `top_n_genes` (or `fdr_threshold` + `effect_size_threshold` when `top_n_genes` is `None`). When `anchor_col` is provided, pseudotime is flipped if it is negatively correlated with that numeric obs column, so `UP`/`DOWN` regulation labels stay stable regardless of which trajectory endpoint was chosen as the origin. Lamian-style visualizations are generated: results summary, heatmap, volcano, per-gene curves, per-sample density/curves, and cluster patterns.

**Source:** `sample_trajectory/trajectory_diff_gene.py:936`

## Signature

```python
def run_trajectory_gam_differential_gene_analysis(
    adata: ad.AnnData,
    pseudotime_source: Union[str, pd.DataFrame, Dict],
    *,
    sample_col: str = "sample",
    celltype_col: Optional[str] = "cell_type",
    batch_col: Optional[Union[str, List[str]]] = None,
    n_features_per_celltype: Optional[int] = 2000,
    columns_to_preserve: Optional[Union[str, List[str]]] = None,
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
    anchor_col: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level AnnData; the per-sample pseudobulk is built internally. |
| `pseudotime_source` | str / DataFrame / dict | — | Path to CSV/TSV, DataFrame, or `{sample_id: pseudotime}` dict. |
| `sample_col` | str | `"sample"` | Sample identifier column in `adata.obs`. |
| `celltype_col` | str, optional | `"cell_type"` | Cell-type column used to build per-cell-type pseudobulk features. |
| `batch_col` | str / list, optional | `None` | Batch column(s) preserved/used when building the pseudobulk. |
| `n_features_per_celltype` | int, optional | `2000` | Number of features kept per cell type when assembling the pseudobulk. |
| `columns_to_preserve` | str / list, optional | `None` | Extra obs columns to carry through onto the pseudobulk samples. |
| `pseudotime_col` | str | `"pseudotime"` | Column within `pseudotime_source` holding the values. |
| `covariate_columns` | list, optional | `None` | Additional fixed-effect covariates (batch, age, sex...). |
| `fdr_threshold` | float | `0.01` | FDR cutoff for pseudoDEG selection. |
| `effect_size_threshold` | float | `1.0` | Minimum effect size (used when `top_n_genes` is `None`). |
| `top_n_genes` | int | `100` | Select at most this many top pseudoDEGs by effect size (set `None` to use the FDR + effect-size rule instead). |
| `num_splines` | int | `5` | GAM basis size. |
| `spline_order` | int | `3` | Spline polynomial order. |
| `output_dir` | str | `"trajectory_diff_gene_results_single"` | Base directory for the results. |
| `visualization_gene_list` | list, optional | `None` | Named genes to always visualize, regardless of rank. |
| `generate_visualizations` | bool | `True` | Turn off to skip the visualization bundle. |
| `group_col` | str, optional | `None` | Optional metadata column for group comparisons in visualizations. |
| `n_clusters` | int | `3` | Number of gene clusters in the heatmap. |
| `top_n_genes_for_curves` | int | `20` | Genes shown in per-gene curve plots. |
| `anchor_col` | str, optional | `None` | Numeric obs column to orient pseudotime; flips it if negatively correlated so `UP`/`DOWN` stay stable. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`pd.DataFrame` — per-gene result table with columns `gene`, `pval`, `dev_exp`, `fdr`, `significant`, `effect_size`, `regulation` (`UP`/`DOWN`), and `pseudoDEG`. An empty DataFrame is returned if no genes are successfully fit.

## Output files

Under `{output_dir}/` (the gene tables are timestamped TSVs):

- `gam_all_genes_<timestamp>.tsv` — full per-gene result table.
- `gam_significant_<timestamp>.tsv` — genes below `fdr_threshold`.
- `gam_pseudoDEGs_<timestamp>.tsv` — selected pseudoDEG table.
- `gam_summary_<timestamp>.txt` — run summary (thresholds, counts, UP/DOWN split).
- `differential_gene_result.txt` — text summary of top genes.

When `generate_visualizations=True`, plots are written under `{output_dir}/visualizations/`:

- `01_results_summary.png`
- `02_tde_heatmap.png`
- `03_xde_heatmap.png` (when `group_col` is set)
- `04_volcano_plot.png`
- `05_gene_curves.png`
- `06_sample_density.png`
- `07_cluster_patterns.png`
- `08_sample_level_curves.png`

## Usage

```python
from sampledisco.sample_trajectory.trajectory_diff_gene import (
    run_trajectory_gam_differential_gene_analysis,
)

results_df = run_trajectory_gam_differential_gene_analysis(
    adata=adata_cell,
    pseudotime_source=pseudotime_dict,
    sample_col="sample",
    celltype_col="cell_type",
    pseudotime_col="pseudotime",
    fdr_threshold=0.05,
    effect_size_threshold=1.0,
    top_n_genes=100,
    num_splines=5,
    spline_order=3,
    anchor_col="sev.level",
    output_dir="/results/rna/trajectoryDEG",
)
```
