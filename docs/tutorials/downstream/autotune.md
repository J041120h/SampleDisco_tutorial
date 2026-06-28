# Autotune (α / block-weight selection)

The sample embedding `uns['X_DR_sample']` is a blend of two views: **composition** blocks (built on `Z_clust`) and an **RMD displacement** block (built on `Z_rmd`). The single knob that controls their relative contribution is the RMD weight **α** — composition vs displacement. `run_autotune` selects that α for you.

`compute_sample_embedding` builds the **multi-resolution** composition blocks (coarse / medium / fine) in one shot, so there is no need to re-cluster at a grid of Leiden resolutions. Set a sensible Leiden resolution once in [`cell_types`](../../api/rna/cell_types_linux.md) (default `leiden_cluster_resolution=0.8`); autotune then searches over α and rebuilds the embedding at the winning weighting.

!!! note "Optional step"
    Autotuning takes longer than the rest of the pipeline because it re-builds the sample embedding for every candidate α. Most users keep the defaults (`rmd_weight=0.60`) and skip it.

## Selecting α with `run_autotune`

Let autotune pick the block weighting:

```python
from sampledisco.parameter_selection.autotune import run_autotune

result = run_autotune(
    adata,
    output_dir="sampledisco_demo_output/rna",
    sample_col="sample",
    celltype_col="cell_type",
    cluster_emb_key="Z_clust",
    rmd_emb_key=None,
    modality_col=None,
    batch_col=None,
    grouping_col="sev.level",
    scoring="auto",
    search="bayesian",
    scope="alpha_only",
    seed=42,
)
```

For multi-omics, pass `modality_col="modality"` and optionally `tune_on_modality="RNA"` to score the search against one modality's labels while still building the final embedding on all units.

See the [`run_autotune` API reference](../../api/shared/run_autotune.md) for the full parameter list.

## Config-driven path

In the config-driven wrapper, autotune is enabled with the per-modality flags `*_autotune_enable` (`rna_autotune_enable` / `atac_autotune_enable` / `multiomics_autotune_enable`). The search is controlled by the companion flags `*_autotune_search`, `*_autotune_scoring`, `*_autotune_alpha_bounds`, and `*_autotune_grouping_col`. The recommended path is simply to set the flag in your YAML and run:

```bash
sampledisco -m complex --config config.yaml
```

## Output

**Writes** → `output_dir/` autotune artifacts: the best parameters, the search trace, and the final sample-level AnnData rebuilt at the winning block weighting (with `uns['X_DR_sample']` set). `run_autotune` returns a dict containing the best params and that AnnData.

## Result

The returned dict carries the selected α and the rebuilt sample AnnData, so the chosen weighting is already baked into `uns['X_DR_sample']`. When run through the wrapper, the selected α is applied automatically before [downstream analysis](../../index.md) proceeds.
