# `find_optimal_cell_resolution_linux`

Two-pass grid search for the Leiden clustering resolution that maximizes the CCA correlation between the sample embedding and a phenotype/trajectory column. The algorithm first sweeps a coarse range (`coarse_start`..`coarse_end` by `coarse_step`), re-running clustering and `calculate_sample_embedding` at every resolution; it then zooms in by `±fine_range` at `fine_step` around the coarse winner. Optionally computes permutation-corrected p-values and writes per-resolution trial artifacts for reproducibility. Works for both RNA and ATAC via the `modality` argument.

**Source:** `parameter_selection/gpu_optimal_resolution.py:33`

## Signature

```python
def find_optimal_cell_resolution_linux(
    adata_cell: AnnData,
    adata_sample: AnnData,
    output_dir: str,
    column: str,
    modality: Literal["rna", "atac"] = "rna",
    trajectory_col: str = "sev.level",
    batch_col: Optional[Union[str, List[str]]] = None,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    cell_embedding_column: str = "X_pca",
    cell_embedding_num_pcs: int = 20,
    n_hvg_features: int = 2000,
    sample_embedding_dimension: int = 10,
    harmony_for_proportion: bool = True,
    preserve_cols_in_sample_embedding: Optional[Union[str, List[str]]] = None,
    n_cca_pcs: int = 10,
    compute_corrected_pvalues: bool = True,
    coarse_start: float = 0.1,
    coarse_end: float = 1.0,
    coarse_step: float = 0.1,
    fine_range: float = 0.02,
    fine_step: float = 0.01,
    verbose: bool = True,
) -> Tuple[float, pd.DataFrame]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata_cell` | AnnData | — | Cell-level AnnData (output of preprocessing). |
| `adata_sample` | AnnData | — | Sample-level AnnData (output of preprocessing). |
| `output_dir` | str | — | Parent output directory; results go under `{output_dir}/{MODALITY}_resolution_optimization_{dr_type}/`. |
| `column` | str | — | DR column to optimize, either `"X_DR_expression"` or `"X_DR_proportion"`. |
| `modality` | `"rna"` or `"atac"` | `"rna"` | Chooses the appropriate HVG / LSI path. |
| `trajectory_col` | str | `"sev.level"` | Phenotype column in `adata_sample.obs` that CCA is tested against. |
| `batch_col` | str or list, optional | `None` | Batch column for sample embedding. |
| `sample_col` | str | `"sample"` | Sample identifier. |
| `celltype_col` | str | `"cell_type"` | Target column where new cluster labels are written. |
| `cell_embedding_column` | str | `"X_pca"` | Cell-level embedding used by the neighbourhood graph. |
| `cell_embedding_num_pcs` | int | `20` | Number of PCs/LSI components. |
| `n_hvg_features` | int | `2000` | HVGs/HVFs per cell type during sample embedding. |
| `sample_embedding_dimension` | int | `10` | DR components for the sample embedding. |
| `harmony_for_proportion` | bool | `True` | Harmony on the proportion embedding. |
| `preserve_cols_in_sample_embedding` | str or list, optional | `None` | Metadata carried into the pseudobulk. |
| `n_cca_pcs` | int | `10` | Number of leading PCs used by CCA. |
| `compute_corrected_pvalues` | bool | `True` | Run permutation-corrected p-values. |
| `coarse_start` / `coarse_end` / `coarse_step` | float | `0.1 / 1.0 / 0.1` | Coarse search grid. |
| `fine_range` / `fine_step` | float | `0.02 / 0.01` | Fine search band around the coarse winner. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`Tuple[float, pd.DataFrame]` — `(best_resolution, trials_df)`. `trials_df` has columns including `resolution`, `cca_score`, `pvalue`, and `corrected_pvalue` (when enabled).

## Output files

Under `{output_dir}/{MODALITY}_resolution_optimization_{expression|proportion}/`:

- `resolutions/{res}/` — per-trial AnnData and diagnostics.
- `summary/optimal.h5ad` — AnnData at the winning resolution.
- `summary/trials.csv` — the returned `trials_df`.

## Usage

```python
from genodistance.parameter_selection import find_optimal_cell_resolution_linux

best_res, trials = find_optimal_cell_resolution_linux(
    adata_cell=adata_cluster,
    adata_sample=adata_sample,
    output_dir="/results/rna",
    column="X_DR_expression",
    modality="rna",
    trajectory_col="sev.level",
    n_cca_pcs=10,
    coarse_start=0.1, coarse_end=1.0, coarse_step=0.1,
    fine_range=0.02,  fine_step=0.01,
)
```
