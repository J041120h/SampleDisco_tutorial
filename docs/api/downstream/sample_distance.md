# `sample_distance`

Unified entry point for pairwise sample-distance computation. Supports three families of methods:

- **Vector metrics** on the single sample embedding (`cosine`, `correlation`, `euclidean`, and any `scipy.spatial.distance.pdist` metric in `VALID_PDIST_METRICS`). Acts on `adata.uns['X_DR_sample']`.
- **EMD (Earth Mover's Distance)** on cell-type proportions paired with cell-type centroids in the supplied cell-level embedding.
- **Distributional distances** (`chi_square`, `jensen_shannon`) on cell-type proportions.

Each call writes a distance matrix CSV and heatmap PDF under `{output_dir}/{method}/`.

**Source:** `sample_distance/sample_distance.py:512`

## Signature

```python
def sample_distance(
    adata: AnnData,
    output_dir: str,
    method: str,
    data_type: str = "ATAC",
    grouping_columns: Optional[List[str]] = None,
    summary_csv_path: Optional[str] = None,
    # EMD-specific parameters
    cell_adata: Optional[AnnData] = None,
    cell_type_column: str = "cell_type",
    sample_column: str = "sample",
    embedding_key: Optional[str] = None,
    n_pcs: int = 20,
    proportions: Optional[pd.DataFrame] = None,
    centroids: Optional[Union[pd.DataFrame, np.ndarray]] = None,
    pseudobulk_adata: Optional[AnnData] = None,
) -> Optional[Dict[str, pd.DataFrame]]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Sample-level AnnData carrying the sample embedding in `.uns['X_DR_sample']`. |
| `output_dir` | str | — | Parent directory; results go under `{output_dir}/{method}/`. |
| `method` | str | — | One of the vector metrics (`"cosine"`, `"correlation"`, `"euclidean"`, ...), `"EMD"`, `"chi_square"`, or `"jensen_shannon"`. |
| `data_type` | str | `"ATAC"` | Modality hint (`"RNA"`, `"ATAC"`, or `"multiomics"`); used to resolve the default cell-level embedding key on the EMD path. |
| `grouping_columns` | list, optional | `None` | `.obs` columns used for grouping annotations on the heatmap. |
| `summary_csv_path` | str, optional | `None` | Optional override for the group-summary CSV path. |
| `cell_adata` | AnnData, optional | `None` | Cell-level AnnData (required for EMD and distributional distances). |
| `cell_type_column` | str | `"cell_type"` | Cell-type column in `cell_adata.obs` (EMD path). |
| `sample_column` | str | `"sample"` | Sample column in `cell_adata.obs` (EMD path). |
| `embedding_key` | str, optional | `None` | Key in `cell_adata.obsm` for cell-type centroid computation (EMD only). When `None`, resolves to `Z_clust` (falling back to `X_pca`/`X_lsi`/`X_glue` per `data_type`). |
| `n_pcs` | int | `20` | Number of PCs to use when computing centroids (EMD only). |
| `proportions` | DataFrame, optional | `None` | Precomputed sample × cell-type proportion matrix (EMD only). |
| `centroids` | DataFrame or ndarray, optional | `None` | Precomputed cell-type centroids (EMD only). |
| `pseudobulk_adata` | AnnData, optional | `None` | Explicit pseudobulk to annotate group metadata (EMD only). |

## Returns

- **Vector metrics:** `Dict[str, DataFrame]` with a single key `"sample_DR"`, whose value is a symmetric `sample × sample` distance matrix.
- **EMD:** `Dict[str, DataFrame]` with key `"EMD"`.
- **`chi_square` / `jensen_shannon`:** `None` — these save their results internally.
- **Unknown method:** `None`.

## Output files

Vector metrics write directly under `{output_dir}/{method}/`:

- `distance_matrix_sample_DR.csv`, `sample_DR_coordinates.csv`
- `sample_distance_sample_DR_heatmap.pdf`
- `distance_statistics_summary_{method}.csv`
- Group-summary CSVs when `grouping_columns` is provided.

EMD and the distributional methods write under their own subfolder (`{output_dir}/{method}/`, with `chi_square` / `jensen_shannon` nested in `{output_dir}/chi_square/` and `{output_dir}/jensen_shannon/`).

## Usage

```python
from sampledisco.sample_distance.sample_distance import sample_distance

for method in ["cosine", "correlation", "euclidean"]:
    sample_distance(
        adata=sample_adata,
        output_dir="/results/rna",
        method=method,
        data_type="RNA",
        grouping_columns=["sev.level"],
    )
```
