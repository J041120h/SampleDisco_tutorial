# `TSCAN`

Unsupervised trajectory inference following the TSCAN paper (Ji & Ji, NAR 2016). Clusters samples with a Gaussian mixture (BIC-selected if `n_clusters=None`), builds a minimum spanning tree on cluster centroids, finds the principal (longest) path through the MST, projects each sample onto its nearest edge, and returns a pseudotime score. Good when you have no supervising phenotype and want structure to emerge from the embedding itself.

**Source:** `sampledisco/sample_trajectory/TSCAN.py:728`

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
| `AnnData_sample` | AnnData | — | Sample-level AnnData with the DR matrix in `.uns`. |
| `column` | str | — | Key for the DR matrix in `.uns`; the single sample embedding key `"X_DR_sample"`. |
| `n_clusters` | int, optional | `None` | Number of sample clusters. When `None`, BIC picks automatically (TSCAN default). |
| `output_dir` | str | `"./"` | Writes to `{output_dir}/TSCAN/`. |
| `grouping_columns` | list, optional | `None` | Metadata columns drawn as overlays on the trajectory plots. |
| `verbose` | bool | `False` | Print progress. |
| `origin` | int, optional | `None` | Cluster index to seed the pseudotime ordering. Must be a principal-path endpoint. When `None`, a random endpoint is chosen. |
| `pseudotime_mode` | str | `"rank"` | `"rank"` or `"distance"`: rank-based (TSCAN default) or projection-distance-based pseudotime. |

## Returns

`Dict` — includes:

- `"pseudotime"` → `{"main_path": {sample_id: float}, "branching_paths": {branch_idx: {sample_id: float}}}`
- `"sample_cluster"` → `{cluster_id: [sample_ids]}`
- `"mst_adjacency"` — adjacency matrix of the cluster MST
- `"main_path"` — ordered list of cluster ids along the principal path
- `"branching_paths"` — side branches off the principal path
- `"pca_data"` — sample coordinates in the DR space used for clustering
- `"graph"` — the MST as a NetworkX graph

`AnnData_sample.obs` is also annotated in place with `tscan_pseudotime_main` and `tscan_cluster`.

## Output files

Under `{output_dir}/TSCAN/`:

- `clusters_by_cluster_{column}.png` — samples colored by inferred cluster.
- `clusters_by_grouping_{column}.png` — samples colored by each entry in `grouping_columns` (only when `grouping_columns` is provided).
- `{column}_pseudotime.csv` — per-sample pseudotime with trajectory/branch/cluster columns.

## Usage

```python
from sampledisco.sample_trajectory.TSCAN import TSCAN

result = TSCAN(
    AnnData_sample=pseudo_adata,
    column="X_DR_sample",
    n_clusters=None,
    output_dir="/results/rna",
    grouping_columns=["sev.level"],
    origin=None,
)
pseudotime = result["pseudotime"]["main_path"]
```
