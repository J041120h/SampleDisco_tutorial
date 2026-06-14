# Sample clustering

`cluster` runs K-means on the single sample embedding `obsm['X_DR_sample']` and returns the legacy `(expr_results, prop_results)` tuple — **both slots now contain the same** `{sample_id: cluster_label}` dict (the two-key `X_DR_expression` / `X_DR_proportion` split is gone). The assignments are the standard input to [`proportion_test`](proportion_test.md) and [`raisinfit`](raisin_dge.md) via their `sample_to_clade` arguments.

!!! note "Legacy arguments"
    `use_expression` and `use_proportion` are still accepted for backward compatibility but are **ignored** — there is only one embedding to cluster.

## Call

```python
from sampledisco.sample_clustering.cluster import cluster

sample_clusters, _ = cluster(
    pseudobulk_adata=pseudo_adata,
    output_dir="/results/rna",
    number_of_clusters=4,
    random_state=0,
)
```

## Output

**Writes** → `/results/rna/sample_cluster/`:

- `kmeans_clusters_sample.csv` — sample ↔ cluster table.
- `kmeans_sample_embedding.png` — 2D scatter of the first two embedding dimensions, colored by cluster.

The pseudobulk AnnData also gets a `cluster_sample_kmeans` column in `.obs`.

## Result

![K-means on sample embedding](../../resource/downstream/kmeans_expression_embedding.png)
<div class="figure-caption">Sample clusters on the first two components of the sample embedding.</div>

See the [API page](../../api/downstream/cluster.md) for the full parameter list.
