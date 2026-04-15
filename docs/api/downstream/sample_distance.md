# `sample_distance`

Unified entry point for pairwise sample-distance computation. Supports three families of methods:

- **Vector metrics** on the dimension-reduced embeddings (`cosine`, `correlation`, `euclidean`, and any `scipy.spatial.distance.pdist` metric). Acts on both `X_DR_expression` and `X_DR_proportion` in `pseudobulk_adata.obsm`.
- **EMD (Earth Mover's Distance)** on cell-type proportions paired with cell-type centroids in the supplied cell-level embedding.
- **Distributional distances** (`chi_square`, `jensen_shannon`) on cell-type proportions.

Each call writes a distance matrix CSV and heatmap PDF under `{output_dir}/Sample_distance/{method}/`.

**Source:** `sample_distance/sample_distance.py:642`

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
    embedding_key: str = "X_pca_harmony",
    n_pcs: int = 20,
    proportions: Optional[pd.DataFrame] = None,
    centroids: Optional[Union[pd.DataFrame, np.ndarray]] = None,
    pseudobulk_adata: Optional[AnnData] = None,
) -> Optional[Dict[str, pd.DataFrame]]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Sample-level (pseudobulk) AnnData with DR in `.obsm`. |
| `output_dir` | str | — | Parent directory; results go under `{output_dir}/Sample_distance/{method}/`. |
| `method` | str | — | One of the vector metrics (`"cosine"`, `"correlation"`, `"euclidean"`, ...), `"EMD"`, `"chi_square"`, or `"jensen_shannon"`. |
| `data_type` | str | `"ATAC"` | Priority hint for picking the best DR key in `.obsm` (`"RNA"` or `"ATAC"`). |
| `grouping_columns` | list, optional | `None` | `.obs` columns used for grouping annotations on the heatmap. |
| `summary_csv_path` | str, optional | `None` | Optional override for the group-summary CSV path. |
| `cell_adata` | AnnData, optional | `None` | Cell-level AnnData (required for EMD and distributional distances). |
| `cell_type_column` | str | `"cell_type"` | Cell-type column in `cell_adata.obs` (EMD path). |
| `sample_column` | str | `"sample"` | Sample column in `cell_adata.obs` (EMD path). |
| `embedding_key` | str | `"X_pca_harmony"` | Key in `cell_adata.obsm` for cell-type centroid computation (EMD only). |
| `n_pcs` | int | `20` | Number of PCs to use when computing centroids (EMD only). |
| `proportions` | DataFrame, optional | `None` | Precomputed sample × cell-type proportion matrix (EMD only). |
| `centroids` | DataFrame or ndarray, optional | `None` | Precomputed cell-type centroids (EMD only). |
| `pseudobulk_adata` | AnnData, optional | `None` | Explicit pseudobulk to annotate group metadata (EMD only). |

## Returns

`Dict[str, DataFrame]` — keyed by the DR used (`"expression"`, `"proportion"`), each value is a symmetric `sample × sample` distance matrix. Returns `None` if the method is not recognized.

## Output files

Under `{output_dir}/Sample_distance/{method}/`:

- `distance_matrix_expression_DR.csv`, `distance_matrix_proportion_DR.csv`
- `expression_DR_heatmap_{method}.pdf`, `proportion_DR_heatmap_{method}.pdf`
- Group-summary CSVs when `grouping_columns` is provided.

## Usage

```python
from genodistance.sample_distance import sample_distance

for method in ["cosine", "correlation", "euclidean"]:
    sample_distance(
        adata=pseudo_adata,
        output_dir="/results/rna",
        method=method,
        data_type="RNA",
        grouping_columns=["sev.level"],
    )
```
