# Autotune (α / block-weight selection)

The sample embedding `uns['X_DR_sample']` is a blend of two views: **composition** blocks (built on `Z_comp`) and an **RMD displacement** block (built on `Z_rmd`). The single knob that controls their relative contribution is the RMD weight **α** — composition vs displacement. `run_autotune` selects that α for you.

`compute_sample_embedding` builds the **multi-resolution** composition blocks (coarse / medium / fine) in one shot, so there is no need to re-cluster at a grid of Leiden resolutions. Set a sensible Leiden resolution once in [`cell_types`](../../api/rna/cell_types_linux.md) (default `leiden_cluster_resolution=0.8`); autotune then searches over α and rebuilds the embedding at the winning weighting.

!!! note "Optional step"
    Autotuning takes longer than the rest of the pipeline because it re-builds the sample embedding for every candidate α. Most users keep the default (`rmd_weight="equal"`: the RMD block gets the same energy as the composition blocks together) and skip it.

## Selecting α with `run_autotune`

Let autotune pick the block weighting:

```python
from sampledisco.parameter_selection.autotune import run_autotune

result = run_autotune(
    adata,
    output_dir="sampledisco_demo_output/rna",
    sample_col="sample",
    celltype_col="cell_type",
    comp_emb_key="Z_comp",
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

!!! note "α is searched on a log scale"

    Every strategy (`bayesian`, `golden_section`, `grid`) explores log10(α), and
    `alpha_bounds` defaults to `(0.1, 100.0)`. The upper end matters: on datasets
    whose signal sits in within-cell-type state rather than in composition, the
    optimum can land above 10 — a stimulation time-course tunes to α≈18, which the
    pre-0.3.0 ceiling of 10.0 truncated. Composition-dominated cohorts tune well
    below 1 and are unaffected. If the reported α comes back equal to either
    bound, widen the range and re-run: the optimum is outside the interval you
    searched.

!!! warning "α is chosen against `grouping_col`"

    The objective scores how well the embedding separates `grouping_col`, so an
    evaluation that uses that same column is circular. Judge a tuned embedding on
    a label that did not take part in the search.

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
