# Trajectory — TSCAN (unsupervised)

When you do not have a supervising phenotype, TSCAN infers a trajectory from the embedding alone: cluster samples with a Gaussian mixture (BIC-selected if `n_clusters=None`), build a minimum spanning tree on cluster centroids, pick the longest path through the tree, and order samples along that path.

## Call

```python
from sampledisco.sample_trajectory.TSCAN import TSCAN

tscan_results = TSCAN(
    AnnData_sample=pseudo_adata,
    column="X_DR_sample",
    n_clusters=None,
    output_dir="/results/rna",
    grouping_columns=["sev.level"],
    origin=None,
    pseudotime_mode="rank",
)
```

TSCAN now runs on the **single** sample embedding `uns['X_DR_sample']` — the two-key
`X_DR_expression` / `X_DR_proportion` split is gone, so there is only one trajectory to compute.
`TSCAN` returns a results `dict` (clusters, MST, principal path, ordering, and pseudotime); it also
writes `tscan_pseudotime_main` and `tscan_cluster` columns back into `AnnData_sample.obs`.

## Output

**Writes** → `/results/rna/TSCAN/` (`{column}` is the embedding key, here `X_DR_sample`):

- `clusters_by_cluster_{column}.png` — points colored by GMM cluster.
- `clusters_by_grouping_{column}.png` — points colored by each entry in `grouping_columns`.
- `{column}_pseudotime.csv` — per-sample pseudotime (with `trajectory_type`, `branch_id`, and `cluster` columns for the main path and any branches).

## Result

![TSCAN clusters on sample embedding](../../resource/downstream/tscan_clusters_by_cluster_expression.png)
![TSCAN colored by severity](../../resource/downstream/tscan_clusters_by_grouping_expression.png)
<div class="figure-caption">Trajectory on the sample embedding, colored either by the inferred cluster or by the phenotype.</div>

See the [API page](../../api/downstream/trajectory_tscan.md) for the full parameter list.
