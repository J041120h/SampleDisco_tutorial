# Resolution selection

Leiden resolution is the most sensitive knob upstream of the sample embedding: too coarse and phenotype signal averages away, too fine and each cell type becomes a noisy isolated cluster. `find_optimal_cell_resolution_linux` automates the choice by doing a two-pass CCA-guided grid search — it sweeps a coarse range first, then zooms in around the winner, re-clustering and re-embedding at every step and ranking each result by the CCA correlation between the sample embedding and a phenotype column. Run it against `X_DR_expression`, `X_DR_proportion`, or both.

!!! note "Optional step"
    This step is advanced and takes longer than the rest of the pipeline because it re-runs clustering + sample embedding for every candidate resolution. Most users pick a reasonable default (e.g. `0.8`) and skip it.

## Call

```python
from genodistance.parameter_selection import find_optimal_cell_resolution_linux

best_res_expr, trials_expr = find_optimal_cell_resolution_linux(
    adata_cell=adata_cluster,
    adata_sample=adata_sample,
    output_dir="/results/rna",
    column="X_DR_expression",
    modality="rna",
    trajectory_col="sev.level",
    n_cca_pcs=10,
    compute_corrected_pvalues=True,
    coarse_start=0.1,
    coarse_end=1.0,
    coarse_step=0.1,
    fine_range=0.02,
    fine_step=0.01,
    verbose=True,
)
```

Run the same call again with `column="X_DR_proportion"` to optimize the proportion embedding.

## Output

**Writes** → `/results/rna/RNA_resolution_optimization_{expression|proportion}/`:

- `resolutions/{res}/` — per-trial AnnData + CCA diagnostics at each resolution tried.
- `summary/optimal.h5ad` — the AnnData at the winning resolution, ready to be used as the new `adata_sample` or pseudobulk.
- `summary/trials.csv` — full per-resolution table (resolution, CCA score, permutation p-value, corrected p-value).

## Result

`(best_resolution, trials_df)` is returned. The `trials_df` lets you inspect how sharply the CCA score peaks and whether the fine search picked up a meaningful improvement over the coarse winner.

Re-run [`cell_types_linux`](../../api/rna/cell_types_linux.md) with `leiden_cluster_resolution=best_res_expr` and then rebuild the sample embedding to lock the optimal resolution into your pipeline state. See the [API page](../../api/shared/find_optimal_cell_resolution_linux.md) for full parameter details.
