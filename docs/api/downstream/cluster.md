# `cluster`

K-means clustering on the sample-level embeddings. Runs once on `X_DR_expression` and once on `X_DR_proportion` (each independently toggleable) and writes both the CSV mapping and a 2D scatter of the first two components colored by cluster. The cluster assignments can be fed directly into [`proportion_test`](proportion_test.md) and [`raisinfit`](raisinfit.md) via the `sample_to_clade` argument.

**Source:** `cluster.py:11`

## Signature

```python
def cluster(
    pseudobulk_adata: ad.AnnData,
    output_dir: str,
    number_of_clusters: int = 5,
    use_expression: bool = True,
    use_proportion: bool = True,
    random_state: int = 0,
) -> Tuple[Optional[Dict[str, int]], Optional[Dict[str, int]]]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudobulk_adata` | AnnData | — | Sample-level AnnData with `X_DR_expression` and/or `X_DR_proportion` in `.obsm`. |
| `output_dir` | str | — | Writes under `{output_dir}/sample_cluster/`. |
| `number_of_clusters` | int | `5` | K for K-means. |
| `use_expression` | bool | `True` | Cluster on `X_DR_expression`. |
| `use_proportion` | bool | `True` | Cluster on `X_DR_proportion`. |
| `random_state` | int | `0` | Seed for reproducibility. |

## Returns

`Tuple[Dict, Dict]` — `(expr_results, prop_results)`, each mapping `{sample_id: cluster_label}`. Either element is `None` if the corresponding flag was disabled.

## Output files

Under `{output_dir}/sample_cluster/`:

- `kmeans_clusters_expression.csv`
- `kmeans_clusters_proportion.csv`
- `kmeans_expression_embedding.png`
- `kmeans_proportion_embedding.png`

The returned AnnData also gets `cluster_expression_kmeans` and `cluster_proportion_kmeans` columns in `.obs`.

## Usage

```python
from genodistance import cluster

expr_clusters, prop_clusters = cluster(
    pseudobulk_adata=pseudo_adata,
    output_dir="/results/rna",
    number_of_clusters=4,
)
```
