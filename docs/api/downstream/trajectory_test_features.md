# `trajectory_test_features`

Runs the trajectory test of [`run_trajectory_gam_differential_gene_analysis`](trajectory_gam_dge.md) on a matrix you have already computed, one row per unit and one column per feature, for example chromVAR motif deviations per sample or cell-type proportions per sample. `adata.X` is used as given: there is no pseudobulk step and no variance-ratio filter (features with variance below 1e-6 are still dropped), and NaN marks a missing unit. Each feature is tested along pseudotime with `test="linear"` (least squares, t-test on the slope; the default here) or `test="gam"` (spline GAM, Wald test of equal spline coefficients), with `batch_col` as an unpenalized fixed effect. BH FDR, effect size, direction (`UP` / `DOWN`), pseudoDEG selection and the output files are those of the trajectory DGE.

Each row needs its own pseudotime, keyed by `adata.obs_names`, and the test treats rows as independent. Run it per cell type (one row per sample), or put cell types in the columns; never test several rows of the same sample together, or the sample is counted several times.

**Import:** `from sampledisco.sample_trajectory.trajectory_diff_gene import trajectory_test_features`

**Source:** `sample_trajectory/trajectory_diff_gene.py:1435`

## Signature

```python
def trajectory_test_features(
    adata: ad.AnnData,
    pseudotime_source: Union[str, pd.DataFrame, Dict],
    *,
    test: str = "linear",
    sample_col: str = "sample",
    batch_col: Optional[Union[str, List[str]]] = None,
    covariate_columns: Optional[List[str]] = None,
    pseudotime_col: str = "pseudotime",
    fdr_threshold: float = 0.05,
    effect_size_threshold: float = 1.0,
    top_n_genes: int = 100,
    num_splines: int = 5,
    spline_order: int = 3,
    n_jobs: int = -1,
    output_dir: str = "trajectory_test_features_results",
    generate_visualizations: bool = False,
    group_col: Optional[str] = None,
    n_clusters: int = 3,
    top_n_genes_for_curves: int = 20,
    anchor_col: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame
```

## Parameters

Parameters shared with [`run_trajectory_gam_differential_gene_analysis`](trajectory_gam_dge.md#parameters) have the same meaning there. The differences:

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | One row per unit (`obs_names` = unit IDs), one column per feature. Used as given. |
| `pseudotime_source` | str, DataFrame or dict | — | Pseudotime keyed by `adata.obs_names`. |
| `test` | str | `"linear"` | `"linear"` or `"gam"`. |
| `batch_col` | str, list or None | `None` | Column(s) of `adata.obs`; an unpenalized fixed effect. |
| `generate_visualizations` | bool | `False` | Write the trajectory-DGE plot set. |

## Returns

`pd.DataFrame`, one row per tested feature, with the columns of the trajectory DGE (`gene` holds the feature name).

## Usage

Cell-type proportions per sample along the CCA pseudotime of the RNA chain, with the batch as a covariate:

```python
import anndata as ad
import pandas as pd
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_trajectory.CCA import CCA_Call
from sampledisco.sample_trajectory.trajectory_diff_gene import trajectory_test_features

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")
_, _, ptime, _ = CCA_Call(pseudo, trajectory_col="sev.level", n_components=2)

prop = pd.crosstab(adata.obs["sample"], adata.obs["cell_type"], normalize="index")
units = ad.AnnData(prop.to_numpy(), obs=pd.DataFrame(index=prop.index.astype(str)),
                   var=pd.DataFrame(index=prop.columns.astype(str)))
units.obs["batch"] = (adata.obs.groupby("sample", observed=True)["batch"].first()
                      .loc[units.obs_names].astype(str).values)

res = trajectory_test_features(units, ptime, test="linear", batch_col="batch",
                               output_dir="api_demo/rna/trajectory_proportions", verbose=False)
print(res[["gene", "pval", "fdr", "effect_size", "regulation"]].sort_values("pval").round(3).to_string(index=False))
```

Output (under a second):

```text
gene  pval   fdr  effect_size regulation
   5 0.009 0.041        5.466       DOWN
   4 0.010 0.041        5.270         UP
   2 0.019 0.051          NaN        NaN
   6 0.253 0.506          NaN        NaN
   8 0.400 0.553          NaN        NaN
   3 0.415 0.553          NaN        NaN
   1 0.549 0.627          NaN        NaN
   7 0.896 0.896          NaN        NaN
```

The demo cell types are Leiden clusters `1` to `8`. Effect size and direction are reported for features with FDR < `fdr_threshold`. Proportions sum to one per sample, so they are not independent features; with 8 samples treat the result as a demonstration of the call, not a finding.
