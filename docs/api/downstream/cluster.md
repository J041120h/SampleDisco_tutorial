# `cluster`

K-means clustering of the samples in the sample embedding. It reads `X_DR_sample` from `.uns` (the `DataFrame`), falling back to `.obsm`, runs scikit-learn `KMeans` on all PCs, adds the labels to `.obs['cluster_sample_kmeans']`, and writes a CSV and a scatter of the first two PCs. The mapping can be passed to [`proportion_test`](proportion_test.md) and [`raisinfit`](raisinfit.md) as `sample_to_clade`.

**Import:** `from sampledisco.sample_clustering.cluster import cluster`

**Source:** `sample_clustering/cluster.py:11`

!!! note "Single-key embedding"
    `use_expression` and `use_proportion` are accepted for backward compatibility and ignored. Both slots of the returned tuple hold the same mapping.

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
| `pseudobulk_adata` | AnnData | — | Sample-level AnnData with `X_DR_sample` in `.uns` or `.obsm`, rows in the same order as the embedding (true for [`build_sample_adata`](../shared/build_sample_adata.md)). Gets `.obs['cluster_sample_kmeans']`. |
| `output_dir` | str | — | Output goes to `{output_dir}/sample_cluster/`. |
| `number_of_clusters` | int | `5` | k. The wrapper's `*_cluster_number` defaults to 4. Keep k well below the number of samples, or clusters of one sample follow. |
| `use_expression`, `use_proportion` | bool | `True` | Ignored. |
| `random_state` | int | `0` | K-means seed. |

## Returns

`(mapping, mapping)`: `{sample_id: cluster}` with integer labels starting at 0, twice.

## Output files

Under `{output_dir}/sample_cluster/`:

- `kmeans_clusters_sample.csv`: columns `sample`, `cluster_sample_kmeans`.
- `kmeans_sample_embedding.png`: PC1 vs PC2 colored by cluster.

## Usage

Continuing from the RNA chain:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_clustering.cluster import cluster

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

clusters, _ = cluster(pseudo, "api_demo/rna", number_of_clusters=2)   # {sample: 0 or 1}
print(pseudo.obs[["sev.level", "cluster_sample_kmeans"]])
```

**Writes** → `api_demo/rna/sample_cluster/kmeans_clusters_sample.csv` and `kmeans_sample_embedding.png` (under a second). The clusters do not follow severity on this demo:

```text
             sev.level cluster_sample_kmeans
sample
CoV-3-Wilk           4                     0
CoV-34-Aruna         4                     1
CoV-36-Aruna         4                     0
CoV-4-Wilk           4                     1
HD-1-Wilk            1                     0
HD-20-Aruna          1                     0
HD-21-Aruna          1                     0
HD-4-Wilk            1                     0
```
