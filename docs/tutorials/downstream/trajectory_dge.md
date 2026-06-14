# Trajectory DGE

Given a pseudotime (from [CCA](trajectory_cca.md) or [TSCAN](trajectory_tscan.md)), this step fits a Generalized Additive Model per gene and ranks by effect size + FDR. It emits Lamian-style visualizations: heatmap, volcano, per-gene curves, sample density, per-sample curves, cluster patterns, and a one-page summary.

## Call

The pseudobulk is built **on the fly** from the cell-level `AnnData` (aggregating cells per sample × cell type), so you pass the preprocessed cell-level `adata` — not a separate sample-level object. `pseudotime_source` accepts a dict, a DataFrame, or a path to a CSV/TSV emitted by [CCA](trajectory_cca.md) or [TSCAN](trajectory_tscan.md).

```python
from sampledisco.sample_trajectory.trajectory_diff_gene import (
    run_trajectory_gam_differential_gene_analysis,
)

results_df = run_trajectory_gam_differential_gene_analysis(
    adata,                               # cell-level preprocessed AnnData
    pseudotime_source=expr_pseudotime,   # dict, DataFrame, or CSV/TSV path from CCA/TSCAN
    sample_col="sample",
    celltype_col="cell_type",
    pseudotime_col="pseudotime",
    covariate_columns=None,
    fdr_threshold=0.01,
    effect_size_threshold=1.0,
    top_n_genes=100,
    num_splines=5,
    spline_order=3,
    output_dir="/results/rna/trajectoryDEG/expression",
    generate_visualizations=True,
    n_clusters=3,
    top_n_genes_for_curves=20,
    anchor_col=None,                     # orient pseudotime so UP/DOWN are stable
)
```

## Output

**Writes** → `/results/rna/trajectoryDEG/expression/`:

| File | Shows |
| --- | --- |
| `pseudoDEGs.csv` | Filtered pseudoDEG table (gene, effect size, p-value, FDR, rank). |
| `all_gene_results.csv` | Full per-gene ranking before filtering. |
| `visualizations/tde_heatmap.png` | Top-gene expression across pseudotime-ordered samples. |
| `visualizations/gene_curves.png` | GAM fits for the top genes. |
| `visualizations/volcano_plot.png` | Effect size vs −log10(FDR). |
| `visualizations/sample_density.png` | Sample distribution along pseudotime. |
| `visualizations/sample_level_curves.png` | Per-sample raw + fitted curves. |
| `visualizations/cluster_patterns.png` | Grouped curves per gene cluster. |
| `visualizations/results_summary.png` | One-page diagnostic. |

## Result

![Top trajectory-differential genes](../../resource/downstream/tde_heatmap.png)
![GAM fits for top genes](../../resource/downstream/gene_curves.png)
![Volcano](../../resource/downstream/volcano_plot.png)
![Sample density along pseudotime](../../resource/downstream/sample_density.png)
![Per-sample trajectories](../../resource/downstream/sample_level_curves.png)
![Cluster-level expression patterns](../../resource/downstream/cluster_patterns.png)
![Results summary](../../resource/downstream/results_summary.png)
<div class="figure-caption">Visualizations from the GAM-based trajectory DGE step.</div>

See the [API page](../../api/downstream/trajectory_gam_dge.md) for the full parameter list.
