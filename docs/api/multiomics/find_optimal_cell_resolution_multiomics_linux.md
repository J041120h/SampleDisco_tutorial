# `find_optimal_cell_resolution_multiomics_linux`

!!! danger "Removed"
    `find_optimal_cell_resolution_multiomics_linux` has been **removed** from `sampledisco`. The CCA-driven cell-resolution sweep for the integrated RNA + ATAC object no longer exists, and there is no drop-in replacement that re-clusters at many Leiden resolutions. Parameter selection is now an **alpha / block-weight autotune** that does not touch the clustering resolution. This page is retained only as a migration guide.

## What replaced it

The old function rebuilt the multi-omics sample embedding at every Leiden resolution and kept the resolution that maximized CCA correlation with a trajectory column. The current pipeline fixes the clustering and instead tunes how the composition and RMD blocks are blended, via [`run_autotune`](#current-approach-run_autotune). For most users this runs automatically inside the config-driven wrapper — you do not call it directly.

### Preferred path — the wrapper

Multi-omics autotune is enabled with a single flag in your config:

```yaml
multiomics_autotune_enable: true
```

Then run the pipeline as usual:

```bash
sampledisco -m complex --config config.yaml
```

The wrapper preprocesses each modality, builds the integrated embedding (`Z_clust` sample-removed, `Z_rmd` sample-preserved), runs the autotune to pick the block-blending weight, and writes the final sample embedding to `adata.uns['X_DR_sample']`.

## Current approach: `run_autotune`

For direct control, call the autotune yourself. It searches the composition-vs-RMD blend (`alpha` / block weights) rather than the clustering resolution.

**Source:** `parameter_selection/autotune.py`

### Signature

```python
from sampledisco.parameter_selection.autotune import run_autotune, DEFAULT_ALPHA_BOUNDS

def run_autotune(
    adata: AnnData,
    output_dir: str,
    *,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    cluster_emb_key: str = "Z_clust",
    rmd_emb_key: Optional[str] = None,
    modality_col: Optional[str] = None,
    batch_col: Optional[Union[str, List[str]]] = None,
    grouping_col: Optional[str] = None,
    medium_K: int = 120,
    fine_K: int = 300,
    rmd_dim: int = 8,
    pca_components: int = 10,
    batch_method: str = "harmony",
    scoring: str = "auto",
    search: str = "bayesian",
    scope: str = "alpha_only",
    alpha_bounds: Tuple[float, float] = DEFAULT_ALPHA_BOUNDS,  # (0.1, 10.0)
    seed: int = 42,
    save: bool = True,
    verbose: bool = True,
    tune_on_modality: Optional[str] = None,
) -> Dict
```

### Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Integrated multi-omics object carrying `Z_clust` and `Z_rmd`. |
| `output_dir` | str | — | Where autotune artifacts are written. |
| `sample_col` | str | `"sample"` | Sample identifier. |
| `celltype_col` | str | `"cell_type"` | Cell-type label column. |
| `cluster_emb_key` | str | `"Z_clust"` | Sample-removed embedding for the composition blocks. |
| `rmd_emb_key` | str, optional | `None` | Sample-preserved (RMD) embedding; resolved automatically when `None`. |
| `modality_col` | str, optional | `None` | Column identifying modality — set to `"modality"` for multi-omics. |
| `batch_col` | str or list, optional | `None` | Batch column(s) for the sample-level correction. |
| `grouping_col` | str, optional | `None` | Grouping used by the leave-one-out RMD reference. |
| `medium_K` / `fine_K` | int | `120 / 300` | Medium- and fine-resolution composition block sizes. |
| `rmd_dim` | int | `8` | RMD displacement dimensions per cluster. |
| `pca_components` | int | `10` | PCA components for the final sample embedding. |
| `batch_method` | str | `"harmony"` | Sample-level batch-correction method. |
| `scoring` | str | `"auto"` | Objective used to score candidate blends. |
| `search` | str | `"bayesian"` | Search strategy over the blend space. |
| `scope` | str | `"alpha_only"` | What to tune (alpha-only by default). |
| `alpha_bounds` | tuple | `(0.1, 10.0)` | Bounds on the block-blend weight. |
| `seed` | int | `42` | Random seed. |
| `save` | bool | `True` | Persist autotune artifacts. |
| `verbose` | bool | `True` | Print progress. |
| `tune_on_modality` | str, optional | `None` | Restrict the bio-preservation/batch proxies to one modality (e.g. `"RNA"`) during the search, while the final embedding is still built on all units. |

### Returns

`Dict` — the best parameters plus the final sample-level AnnData (with `uns['X_DR_sample']` populated).

## Usage

```python
from sampledisco.parameter_selection.autotune import run_autotune

result = run_autotune(
    adata_integrated,
    output_dir="/results/multiomics",
    modality_col="modality",
    tune_on_modality="RNA",
)
```

To fold a previously selected optimal embedding back into the merged object, see [`replace_optimal_dimension_reduction`](replace_optimal_dimension_reduction.md).
