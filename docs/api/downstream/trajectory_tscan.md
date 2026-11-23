# `TSCAN`

Unsupervised sample trajectory, following TSCAN (Ji & Ji, *Nucleic Acids Res.* 2016):

1. Cluster the samples in the embedding with a Gaussian mixture (full covariance). With `n_clusters=None` the number is chosen by BIC over k = 2 … min(20, n_samples // 2).
2. Build a minimum spanning tree (MST) on the cluster centroids and take its longest path as the **main path**; the other branches become branching paths.
3. Project each sample onto the nearest MST edge and order the samples along each path.
4. Pseudotime is either the rank along the path scaled to [0, 1] (`"rank"`, the TSCAN default) or the cumulative Euclidean distance along the path scaled to [0, 1] (`"distance"`).

The path starts at `origin`. With `origin=None` the start is chosen deterministically: the endpoint cluster with the lower mean of the first **numeric** column in `grouping_columns`, or else the endpoint with the lower index.

**Import:** `from sampledisco.sample_trajectory.TSCAN import TSCAN`

**Source:** `sample_trajectory/TSCAN.py:728`

## Signature

```python
def TSCAN(
    AnnData_sample: sc.AnnData,
    column: str,
    n_clusters: Optional[int] = None,
    output_dir: str = "./",
    grouping_columns: Optional[List[str]] = None,
    verbose: bool = False,
    origin: Optional[int] = None,
    pseudotime_mode: str = "rank",
) -> Dict
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `AnnData_sample` | AnnData | — | Sample-level AnnData; `.uns[column]` must be a `DataFrame` (as from [`build_sample_adata`](../shared/build_sample_adata.md)). |
| `column` | str | — | `.uns` key of the embedding, normally `"X_DR_sample"`. |
| `n_clusters` | int or None | `None` | Number of sample clusters; `None` = BIC. Wrapper key: `*_tscan_n_clusters`. |
| `output_dir` | str | `"./"` | Output goes to `{output_dir}/TSCAN/`. |
| `grouping_columns` | list or None | `None` | `.obs` columns drawn on the plots; the first numeric one also orients the path when `origin=None`. Wrapper key: `*_trajectory_visualization_label`. |
| `verbose` | bool | `False` | Print progress. |
| `origin` | int or None | `None` | Index (0-based) of the endpoint cluster to start from; it must be an end of the main path, otherwise `ValueError`. Wrapper key: `*_tscan_origin`. |
| `pseudotime_mode` | str | `"rank"` | `"rank"` or `"distance"`. Wrapper key: `*_tscan_pseudotime_mode`. |

## Returns

`dict`, including:

| Key | Value |
| --- | --- |
| `"pseudotime"` | `{"main_path": {sample: t}, "branching_paths": {branch_index: {sample: t}}}` |
| `"main_path"`, `"branching_paths"` | Cluster indices along the main path and each branch. |
| `"origin"` | The starting cluster index. |
| `"sample_cluster"` | `{cluster_name: [samples]}` |
| `"ordered_samples"`, `"sample_projections"` | Main-path order and per-sample projections (and `*_branching_paths` versions). |
| `"mst_adjacency"`, `"graph"` | The MST as a matrix and as a NetworkX graph. |
| `"pca_data"` | The embedding used. |
| `"cluster_names"`, `"n_samples_total"` | Bookkeeping. |

`AnnData_sample.obs` also gets `tscan_pseudotime_main` and `tscan_cluster`.

## Output files

Under `{output_dir}/TSCAN/`:

- `clusters_by_cluster_{column}.png`: samples colored by cluster, with the MST.
- `clusters_by_grouping_{column}.png`: samples colored by the grouping columns (only with `grouping_columns`).
- `{column}_pseudotime.csv`: `sample_id`, `trajectory_type`, `branch_id`, `pseudotime`, `cluster`.

## Usage

Continuing from the RNA chain:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_trajectory.TSCAN import TSCAN

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

res = TSCAN(pseudo, "X_DR_sample", output_dir="api_demo/rna", grouping_columns=["sev.level"])
print(res["main_path"], res["origin"])
print(pseudo.obs[["sev.level", "tscan_cluster", "tscan_pseudotime_main"]].sort_values("tscan_pseudotime_main"))
```

**Writes** → `api_demo/rna/TSCAN/clusters_by_cluster_X_DR_sample.png`, `clusters_by_grouping_X_DR_sample.png` and `X_DR_sample_pseudotime.csv`. BIC picks 2 clusters, so the main path is a single edge (about 4 s):

```text
[0, 1] 0
             sev.level tscan_cluster  tscan_pseudotime_main
sample
HD-4-Wilk            1     cluster_1               0.000000
HD-1-Wilk            1     cluster_1               0.142857
HD-21-Aruna          1     cluster_1               0.285714
HD-20-Aruna          1     cluster_1               0.428571
CoV-36-Aruna         4     cluster_2               0.571429
CoV-3-Wilk           4     cluster_2               0.714286
CoV-4-Wilk           4     cluster_2               0.857143
CoV-34-Aruna         4     cluster_2               1.000000
```

The two clusters match severity here, and the rank pseudotime runs from the healthy cluster to the severe one (in an earlier run on another node, `CoV-3-Wilk` fell in the healthy cluster).
