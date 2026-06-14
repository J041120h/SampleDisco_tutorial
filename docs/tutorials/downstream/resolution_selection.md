# Resolution selection

!!! warning "Removed in the current release"
    The CCA-guided resolution sweep `find_optimal_cell_resolution_linux` (and its multi-omics twin) **no longer exists** in `sampledisco`. The two-key `X_DR_expression` / `X_DR_proportion` embedding it ranked against has been replaced by the single sample embedding `uns['X_DR_sample']`, and parameter selection is now done by **alpha / block-weight autotune** rather than by re-clustering at every candidate resolution. This page is kept as a redirect — see the current approach below.

Leiden resolution still matters upstream of the sample embedding: too coarse and phenotype signal averages away, too fine and each cell type becomes a noisy isolated cluster. But the package no longer sweeps it with a CCA-ranked grid search. Instead, `compute_sample_embedding` builds **multi-resolution** composition blocks (coarse / medium / fine) in one shot, and the tunable knob is the relative weight between the composition and RMD (displacement) blocks. That weight is what `run_autotune` selects.

!!! note "Optional step"
    Autotuning is advanced and takes longer than the rest of the pipeline because it re-builds the sample embedding for every candidate parameter set. Most users keep the defaults (Leiden resolution `0.8`, `rmd_weight=0.60`) and skip it.

## Current approach — autotune

Set a sensible Leiden resolution in [`cell_types`](../../api/rna/cell_types_linux.md) (default `leiden_cluster_resolution=0.8`), then let autotune pick the block weighting:

```python
from sampledisco.parameter_selection.autotune import run_autotune

result = run_autotune(
    adata,
    output_dir="/results/rna",
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

In the config-driven wrapper this is enabled with the per-modality flags `rna_autotune_enable` / `atac_autotune_enable` / `multiomics_autotune_enable`, so the recommended path is simply to set the flag in your YAML and run:

```bash
sampledisco -m complex --config config.yaml
```

## Output

**Writes** → `output_dir/` autotune artifacts: the best parameters, the search trace, and the final sample-level AnnData rebuilt at the winning block weighting (with `uns['X_DR_sample']` set). `run_autotune` returns a dict containing the best params and that AnnData.

## Result

The returned dict carries the selected weighting and the rebuilt sample AnnData, so the optimal parameters are already baked into `uns['X_DR_sample']` — no manual re-clustering loop is needed. When run through the wrapper, the chosen weighting is applied automatically before [downstream analysis](../../index.md) proceeds.
