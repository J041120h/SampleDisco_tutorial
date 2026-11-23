# `run_trajectory_gam_differential_gene_analysis`

Finds genes whose expression changes along a sample pseudotime. It takes a **cell-level** AnnData (normalized, log-transformed `.X`), builds a per-sample pseudobulk internally (per cell type when `celltype_col` is set: features are cell type × gene, keeping up to `n_features_per_celltype` per cell type, with optional batch handling), then fits one generalized additive model (GAM) per feature with pseudotime as a smooth term and `covariate_columns` as fixed effects. Effect sizes and deviance-based p-values are computed per feature, BH-FDR corrected, and **pseudoDEGs** are selected: the `top_n_genes` largest effects among FDR-significant features, or, with `top_n_genes=None`, all features passing `fdr_threshold` and `effect_size_threshold`. With `anchor_col`, the pseudotime is flipped when it correlates negatively with that numeric column, so `UP` / `DOWN` keep a fixed meaning. Lamian-style plots follow.

The spline basis is shrunk automatically for small cohorts, and a feature is fitted only if it is non-zero in at least `spline_order + num_splines + 2` samples (after that adjustment). With 8 samples and the default `spline_order=3` no feature qualifies; use `spline_order=2`, as in the example and the demo config.

**Import:** `from sampledisco.sample_trajectory.trajectory_diff_gene import run_trajectory_gam_differential_gene_analysis`

**Source:** `sample_trajectory/trajectory_diff_gene.py:973`

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
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 1.0,
    top_n_genes: int = 100,
    num_splines: int = 5,
    spline_order: int = 3,
    n_jobs: int = -1,
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
| `adata` | AnnData | — | Cell-level AnnData; the pseudobulk is built internally. For multi-omics use the RNA file (`preprocess/adata_rna_preprocessed.h5ad`), not the embedding-only union. |
| `pseudotime_source` | str, DataFrame or dict | — | `{sample_id: pseudotime}` (e.g. from [`CCA_Call`](trajectory_cca_call.md) or [`TSCAN`](trajectory_tscan.md)), a DataFrame, or a CSV/TSV path. |
| `sample_col` | str | `"sample"` | Sample column. |
| `celltype_col` | str or None | `"cell_type"` | Split the pseudobulk by cell type; `None` averages all cells per sample. |
| `batch_col` | str, list or None | `None` | Batch column(s) used when building the pseudobulk (several are combined into one). |
| `n_features_per_celltype` | int or None | `2000` | Features kept per cell type. |
| `columns_to_preserve` | str, list or None | `None` | Extra `.obs` columns carried onto the pseudobulk samples. |
| `pseudotime_col` | str | `"pseudotime"` | Pseudotime column when the source is a table. |
| `covariate_columns` | list or None | `None` | Fixed-effect covariates (sample-level `.obs` columns). Wrapper key: `*_trajectory_diff_gene_covariate`. |
| `fdr_threshold` | float | `0.05` | FDR cutoff. |
| `effect_size_threshold` | float | `1.0` | Minimum effect size (used when `top_n_genes` is `None`). |
| `top_n_genes` | int or None | `100` | Number of pseudoDEGs; `None` switches to the threshold rule. |
| `num_splines` / `spline_order` | int | `5` / `3` | Requested GAM basis; reduced automatically for small cohorts (see above). |
| `n_jobs` | int | `-1` | Processes for the per-feature fits (`-1` = all cores; results are identical to serial). |
| `output_dir` | str | `"trajectory_diff_gene_results_single"` | Output folder. The wrapper uses `{modality dir}/trajectoryDEG`. |
| `visualization_gene_list` | list or None | `None` | Genes that always get their own plots (under `gene_visualizations/`). |
| `generate_visualizations` | bool | `True` | Write the plot set below. |
| `group_col` | str or None | `None` | Group column for the XDE heatmap. |
| `n_clusters` | int | `3` | Gene clusters in the heatmaps. |
| `top_n_genes_for_curves` | int | `20` | Genes shown in the curve plots. |
| `anchor_col` | str or None | `None` | Numeric `.obs` column to orient the pseudotime. The wrapper passes `trajectory_col` with TSCAN and nothing with CCA. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`pd.DataFrame`, one row per fitted feature, with `gene`, `pval`, `dev_exp`, `fdr`, `significant`, `effect_size`, `regulation` (`UP` / `DOWN`) and `pseudoDEG`. Empty if no feature could be fitted.

## Output files

Under `{output_dir}/`, with a timestamp in the table names:

- `gam_all_genes_<timestamp>.tsv`, `gam_significant_<timestamp>.tsv` (FDR < `fdr_threshold`; may be empty), `gam_pseudoDEGs_<timestamp>.tsv` (may be empty), `gam_summary_<timestamp>.txt`.
- `differential_gene_result.txt`: text summary.
- With `generate_visualizations=True`, under `visualizations/`: `01_results_summary.png`, `02_tde_heatmap.png`, `03_xde_heatmap.png` (with `group_col`), `04_volcano_plot.png`, `05_gene_curves.png`, `06_sample_density.png`, `07_cluster_patterns.png`, `08_sample_level_curves.png`. Plots that need significant genes are skipped when there are none.

## Usage

Continuing from the RNA chain, with the CCA pseudotime:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_trajectory.CCA import CCA_Call
from sampledisco.sample_trajectory.trajectory_diff_gene import (
    run_trajectory_gam_differential_gene_analysis,
)

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")
_, _, ptime, _ = CCA_Call(pseudo, trajectory_col="sev.level", n_components=2)

res = run_trajectory_gam_differential_gene_analysis(
    adata,
    ptime,
    sample_col="sample",
    celltype_col="cell_type",
    spline_order=2,        # 8 samples: spline_order=3 leaves nothing to fit
    n_jobs=8,
    output_dir="api_demo/rna/trajectoryDEG",
)
print(len(res), "features fitted;", int(res["significant"].sum()), "with FDR < 0.05")
```

**Writes** → under `api_demo/rna/trajectoryDEG/`: `gam_all_genes_<timestamp>.tsv`, `gam_significant_<timestamp>.tsv`, `gam_pseudoDEGs_<timestamp>.tsv`, `gam_summary_<timestamp>.txt`, `differential_gene_result.txt`, and `visualizations/01_results_summary.png`, `04_volcano_plot.png`, `05_gene_curves.png`, `06_sample_density.png`, `08_sample_level_curves.png` (the heatmaps need significant genes and are skipped). From the log:

```text
[CCA] Flipped pseudotime so it increases with 'sev.level' (Spearman rho before flip: -0.873)
...
[pseudobulk] 8 samples × 8 cell types
[pseudobulk] result shape: (8, 16000)
  → Spline parameters: n_splines=3, order=2
  → Successfully fitted 10033/16000 genes
  → Skipped 5967 genes:
      too_few_samples: 5967
  → Found 0 significant genes (FDR < 0.05)
10033 features fitted; 0 with FDR < 0.05
```

The function also warns (`RuntimeWarning`) that a spline GAM is underpowered below about 12 samples: on 8 samples an empty result says little. About 18 s, 3.2 GB including the worker processes.

!!! note "When nothing can be fitted (0.6.0)"
    If no feature passes the identifiability check (for example `spline_order=3` on 8 samples), 0.6.0 raises a `RuntimeWarning` that names the number of non-zero samples needed and the adjusted spline settings, and still returns a DataFrame (empty) and writes the tables. It does not return `None`.
