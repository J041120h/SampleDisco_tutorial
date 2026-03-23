# `run_trajectory_gam_differential_gene_analysis`

**Module**: `code/sample_trajectory/trajectory_diff_gene.py`

```python
run_trajectory_gam_differential_gene_analysis(
    pseudobulk_adata,
    pseudotime_source,
    *,
    sample_col="sample",
    pseudotime_col="pseudotime",
    covariate_columns=None,
    fdr_threshold=0.01,
    effect_size_threshold=1.0,
    top_n_genes=100,
    num_splines=5,
    spline_order=3,
    output_dir="trajectory_diff_gene_results_single",
    visualization_gene_list=None,
    generate_visualizations=True,
    group_col=None,
    n_clusters=3,
    top_n_genes_for_curves=20,
    verbose=True,
)
```

GAM-based trajectory differential analysis for one pseudotime vector.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `pseudobulk_adata` | required | Pseudobulk sample x gene AnnData. |
| `pseudotime_source` | required | Pseudotime input (file/DataFrame/dict). |
| `sample_col` | `"sample"` | Sample ID column. |
| `pseudotime_col` | `"pseudotime"` | Pseudotime value column. |
| `covariate_columns` | `None` | Additional covariates in GAM model. |
| `fdr_threshold` | `0.01` | FDR cutoff for significance. |
| `effect_size_threshold` | `1.0` | Effect-size threshold for pseudoDEG flag. |
| `top_n_genes` | `100` | Maximum pseudoDEGs retained. |
| `num_splines` | `5` | Number of spline basis functions. |
| `spline_order` | `3` | Spline polynomial order. |
| `output_dir` | `"trajectory_diff_gene_results_single"` | Output directory. |
| `visualization_gene_list` | `None` | Optional explicit genes for curve plotting. |
| `generate_visualizations` | `True` | Generate result summary plots. |
| `group_col` | `None` | Group column for grouped visualizations. |
| `n_clusters` | `3` | Cluster count used in heatmap patterning. |
| `top_n_genes_for_curves` | `20` | Number of genes used in curve panels. |
| `verbose` | `True` | Print progress logs. |

## Returns

- DataFrame of gene statistics and pseudoDEG labels.

