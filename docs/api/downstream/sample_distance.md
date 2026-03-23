# `sample_distance`

**Module**: `code/sample_distance/sample_distance.py`

```python
sample_distance(
    adata,
    output_dir,
    method,
    data_type="ATAC",
    grouping_columns=None,
    summary_csv_path=None,
    cell_adata=None,
    cell_type_column="cell_type",
    sample_column="sample",
    embedding_key="X_pca_harmony",
    n_pcs=20,
    proportions=None,
    centroids=None,
    pseudobulk_adata=None,
)
```

Unified sample distance API across DR metrics and composition-based methods.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Sample-level AnnData (or metadata source). |
| `output_dir` | required | Output path for distance results. |
| `method` | required | Distance method (`cosine`, `correlation`, `EMD`, `chi_square`, `jensen_shannon`, etc.). |
| `data_type` | `"ATAC"` | DR key prioritization behavior (`ATAC`/`RNA`). |
| `grouping_columns` | `None` | Group columns used for summary/statistics. |
| `summary_csv_path` | `None` | CSV path for summary output. |
| `cell_adata` | `None` | Cell-level AnnData (required for EMD/chi-square/JS). |
| `cell_type_column` | `"cell_type"` | Cell type label column. |
| `sample_column` | `"sample"` | Sample ID column. |
| `embedding_key` | `"X_pca_harmony"` | Embedding key used in EMD mode. |
| `n_pcs` | `20` | Number of embedding PCs used by EMD. |
| `proportions` | `None` | Optional precomputed cell-type proportions for EMD. |
| `centroids` | `None` | Optional precomputed cell-type centroids for EMD. |
| `pseudobulk_adata` | `None` | Pseudobulk AnnData for metadata lookup. |

