# `run_autotune`

Selects the **RMD-weight α** that best balances the cellular-composition blocks against the reference-relative mean-displacement (RMD) block when building the sample embedding. This is SampleDisco's parameter-selection step — it replaces the removed CCA-guided cell-resolution sweep.

Given a preprocessed cell-level `AnnData` (carrying `Z_clust` and `Z_rmd`), autotune builds the composition + RMD blocks **once**, then searches over α (grid / golden-section / Bayesian) to maximize an adaptive proxy objective — phenotype alignment (CCA / PC-R²) combined with batch-mixing scores — and rebuilds `uns['X_DR_sample']` at the winning α.

**Source:** `parameter_selection/autotune.py:543`

```python
from sampledisco.parameter_selection.autotune import run_autotune

result = run_autotune(
    adata,
    output_dir,
    sample_col="sample",
    celltype_col="cell_type",
    cluster_emb_key="Z_clust",
    grouping_col="sev.level",   # phenotype to align alpha against
    scope="alpha_only",
    search="bayesian",          # "bayesian" | "golden" | "grid"
    scoring="auto",
    alpha_bounds=(0.1, 10.0),
    seed=42,
)
```

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Preprocessed cell-level object carrying `Z_clust` and `Z_rmd`. |
| `output_dir` | str | — | Where autotune artifacts (best params, search trace, rebuilt AnnData) are written. |
| `grouping_col` | str | `None` | Phenotype/condition column the search aligns α to (CCA / PC-R²). |
| `cluster_emb_key` / `rmd_emb_key` | str | `"Z_clust"` / `None` | Composition (sample-removed) and RMD (sample-preserved) embedding keys. |
| `scope` | str | `"alpha_only"` | Search scope. Only `"alpha_only"` is currently supported. |
| `search` | str | `"bayesian"` | Strategy: `bayesian` (GP), `golden` (golden-section), or `grid`. |
| `scoring` | str | `"auto"` | Proxy objective; `auto` selects an ensemble from the available metadata. |
| `alpha_bounds` | tuple | `(0.1, 10.0)` | Lower / upper bounds for α. |
| `tune_on_modality` | str | `None` | (Multi-omics) restrict the scoring proxies to one modality's units while still building the final embedding on **all** units. |
| `medium_K` / `fine_K` / `rmd_dim` / `pca_components` | int | `120` / `300` / `8` / `10` | Block-construction knobs (mirror `compute_sample_embedding`). |
| `batch_method` | str | `"harmony"` | Sample-level batch correction used when rebuilding the embedding. |
| `seed` | int | `42` | RNG seed. |
| `save` | bool | `True` | Write artifacts + rebuilt embedding to `output_dir`. |

**Returns** — a `dict` containing the best parameters, the search trace, and the final sample-level `AnnData` rebuilt at the winning α (with `uns['X_DR_sample']` set).

!!! tip "Use the config-driven wrapper"
    Most users enable autotune through the per-modality flags `rna_autotune_enable` / `atac_autotune_enable` / `multiomics_autotune_enable` (with `*_autotune_search`, `*_autotune_scoring`, `*_autotune_alpha_bounds`, `*_autotune_grouping_col`) in the YAML config, then run `sampledisco -m complex --config config.yaml`. The wrapper calls `run_autotune` in place of `compute_sample_embedding` when the flag is set. See the [autotune walkthrough](../../tutorials/downstream/resolution_selection.md).
