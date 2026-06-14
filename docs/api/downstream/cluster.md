# `cluster`

K-means clustering on the single-key sample-level embedding `obsm['X_DR_sample']`. Writes the CSV mapping of samples to clusters and a 2D scatter of the first two embedding dimensions colored by cluster. The cluster assignments can be fed directly into [`proportion_test`](proportion_test.md) and [`raisinfit`](raisinfit.md) via the `sample_to_clade` argument.

!!! note "Single-key embedding"
    The current pipeline produces **one** sample embedding, `X_DR_sample`. The legacy `use_expression` / `use_proportion` flags and the old `X_DR_expression` / `X_DR_proportion` keys are gone — both flags are accepted for backward compatibility but **ignored**, and the returned tuple has the same mapping in both slots.

**Source:** `sampledisco/sample_clustering/cluster.py:11`

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
| `pseudobulk_adata` | AnnData | — | Sample-level AnnData with `X_DR_sample` in `.obsm`. |
| `output_dir` | str | — | Writes under `{output_dir}/sample_cluster/`. |
| `number_of_clusters` | int | `5` | K for K-means. |
| `use_expression` | bool | `True` | Legacy argument — **accepted but ignored**. |
| `use_proportion` | bool | `True` | Legacy argument — **accepted but ignored**. |
| `random_state` | int | `0` | Seed for reproducibility. |

## Returns

`Tuple[Dict, Dict]` — `(expr_results, prop_results)`, kept for backward compatibility. Both slots contain the **same** `{sample_id: cluster_label}` mapping computed on `X_DR_sample`.

## Output files

Under `{output_dir}/sample_cluster/`:

- `kmeans_clusters_sample.csv`
- `kmeans_sample_embedding.png`

The passed-in AnnData also gets a `cluster_sample_kmeans` column in `.obs`.

## Usage

```python
from sampledisco.sample_clustering.cluster import cluster

sample_clusters, _ = cluster(
    pseudobulk_adata=pseudo_adata,
    output_dir="/results/rna",
    number_of_clusters=4,
)
```
