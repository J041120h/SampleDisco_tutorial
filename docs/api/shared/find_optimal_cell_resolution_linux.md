# `find_optimal_cell_resolution_linux`

!!! danger "Removed"
    `find_optimal_cell_resolution_linux` has been **removed** from `sampledisco` and has **no drop-in replacement**. The CCA-driven, two-pass Leiden-resolution sweep no longer exists, and the legacy two-key embedding it optimized (`X_DR_expression` / `X_DR_proportion`) is gone — the current pipeline produces a single sample embedding at `adata.uns['X_DR_sample']`. Parameter selection is now an **alpha / block-weight autotune** (`sampledisco.parameter_selection.autotune.run_autotune`), wired into the config-driven wrapper. This page is retained only as a migration pointer.

## What it used to do

The old function ran a two-pass grid search for the Leiden clustering resolution that maximized the CCA correlation between the sample embedding and a phenotype/trajectory column, re-running clustering and the (now-removed) `calculate_sample_embedding` at every resolution. That entire mechanism — coarse/fine resolution sweeps, per-resolution trial artifacts, and the `X_DR_expression`/`X_DR_proportion` columns it tuned — has been deleted.

## What to use instead

Instead of sweeping Leiden resolution, the current package tunes the blend weight (alpha) between the composition blocks and the RMD displacement block of the single-key embedding. Two ways to use it:

### 1. Via the config-driven wrapper (recommended)

Enable autotune per modality with the `*_autotune_enable` flags in your YAML config, then run the standard entry point:

```bash
sampledisco -m complex --config config.yaml
```

```yaml
# config.yaml (excerpt)
rna_autotune_enable: true        # RNA
atac_autotune_enable: true       # ATAC
multiomics_autotune_enable: true # multi-omics
```

The wrapper calls `run_autotune` internally after the sample embedding step and writes the tuned result in place.

### 2. Calling `run_autotune` directly

**Source:** `parameter_selection/autotune.py`

```python
from sampledisco.parameter_selection.autotune import run_autotune

result = run_autotune(
    adata,                       # cell-level AnnData carrying Z_clust / Z_rmd
    output_dir="/results/rna",
    sample_col="sample",
    celltype_col="cell_type",
    cluster_emb_key="Z_clust",
    scoring="auto",
    search="bayesian",
    scope="alpha_only",
)
```

For multi-omics, restrict the search to one modality:

```python
result = run_autotune(
    adata,
    output_dir="/results/multiomics",
    modality_col="modality",
    tune_on_modality="RNA",
)
```

## Signature (`run_autotune`)

```python
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

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level AnnData carrying `Z_clust` (sample-removed) and `Z_rmd` (sample-preserved). |
| `output_dir` | str | — | Output directory for the autotune report and tuned embedding. |
| `sample_col` | str | `"sample"` | Sample identifier column in `adata.obs`. |
| `celltype_col` | str | `"cell_type"` | Cell-type label column. |
| `cluster_emb_key` | str | `"Z_clust"` | Sample-removed cell embedding used for the composition blocks. |
| `rmd_emb_key` | str, optional | `None` | RMD (sample-preserved) embedding key; resolves to `Z_rmd` when `None`. |
| `modality_col` | str, optional | `None` | Modality column (set for multi-omics, e.g. `"modality"`). |
| `batch_col` | str or list, optional | `None` | Batch column(s) for sample-level Harmony correction. |
| `grouping_col` | str, optional | `None` | Optional grouping used by some scoring proxies. |
| `medium_K` / `fine_K` | int | `120` / `300` | Medium- and fine-resolution composition block sizes. |
| `rmd_dim` | int | `8` | RMD displacement dimensions per cluster. |
| `pca_components` | int | `10` | Sample-embedding DR components. |
| `batch_method` | str | `"harmony"` | Sample-level batch-correction method. |
| `scoring` | str | `"auto"` | Scoring objective driving the search. |
| `search` | str | `"bayesian"` | Search strategy (`bayesian`, `grid`, `golden`). |
| `scope` | str | `"alpha_only"` | Tuning scope (currently `alpha_only`). |
| `alpha_bounds` | tuple | `(0.1, 10.0)` | Search bounds for the blend weight alpha. |
| `seed` | int | `42` | Random seed. |
| `save` | bool | `True` | Write the tuned embedding and report. |
| `verbose` | bool | `True` | Print progress. |
| `tune_on_modality` | str, optional | `None` | Restrict the search proxies to one modality (e.g. `"RNA"`); final embedding still built on all units. |

## Returns

`Dict` — the best parameters (including the tuned alpha) plus the final autotuned sample-level result. The tuned embedding is written in place to `adata.uns['X_DR_sample']` (units × PCs).
