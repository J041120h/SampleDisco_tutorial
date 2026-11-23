# `sample_distance`

Pairwise distances between samples. Three families of `method`:

- **Vector metrics**, on the sample embedding `adata.uns['X_DR_sample']`: `euclidean`, `sqeuclidean`, `minkowski`, `cityblock`, `chebyshev`, `cosine`, `correlation`, `hamming`, `jaccard`, `canberra`, `braycurtis`, `matching` (passed to `scipy.spatial.distance.pdist`). The matrix is divided by its maximum, so values lie in [0, 1].
- **`"EMD"`**: Earth Mover's Distance between samples' cell-type proportions, with ground distances between cell-type centroids in a cell-level embedding. Needs `cell_adata`.
- **`"chi_square"`**, **`"jensen_shannon"`**: distances between cell-type proportion vectors. Need `cell_adata`. Broken on 0.5.0; see the warning under Output files.

For every method, `grouping_columns` drives a **group check**: the mean between-group distance divided by the mean within-group distance (above 1 means the groups separate), written to a text file and, with `summary_csv_path`, into a summary CSV. If the grouping cannot be resolved, the check is skipped with a warning and the distance matrix is still written.

**Import:** `from sampledisco.sample_distance.sample_distance import sample_distance`

**Source:** `sample_distance/sample_distance.py:540`

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
| `adata` | AnnData | — | AnnData with `uns['X_DR_sample']`, normally the sample-level one from [`build_sample_adata`](../shared/build_sample_adata.md). Its `.obs` supplies the groups for vector metrics. In 0.6.0 the group check also accepts the cell-level object: it then takes the first grouping value per `obs['sample']` (0.5.0 skipped the check with a warning in that case). |
| `output_dir` | str | — | Parent folder (see Output files). |
| `method` | str | — | A vector metric above, `"EMD"`, `"chi_square"` or `"jensen_shannon"`. Anything else prints a warning and returns `None`. |
| `data_type` | str | `"ATAC"` | `"RNA"`, `"ATAC"` or `"multiomics"`. Only used by EMD to pick the fallback cell embedding when `Z_comp` is missing (`X_pca`, `X_lsi`, `X_glue`). |
| `grouping_columns` | list or None | `None` | `.obs` columns for the group check (the first valid one is used). |
| `summary_csv_path` | str or None | `None` | CSV updated with the group-check score (row = `sample_DR`/`EMD`/…, column = method). |
| `cell_adata` | AnnData or None | `None` | Cell-level AnnData; required for EMD, chi-square and Jensen–Shannon. |
| `cell_type_column` / `sample_column` | str | `"cell_type"` / `"sample"` | Columns in `cell_adata.obs` (EMD, chi-square, Jensen–Shannon). |
| `embedding_key` | str or None | `None` | EMD only: `.obsm` key for the centroids. `None` resolves `Z_comp` (or a legacy name), then the `data_type` fallback. |
| `n_pcs` | int | `20` | EMD only: columns of the embedding used for centroids. |
| `proportions` / `centroids` | DataFrame / array or None | `None` | EMD only: precomputed sample × cell-type proportions and cell-type centroids. |
| `pseudobulk_adata` | AnnData or None | `None` | Sample-level AnnData whose `.obs` supplies the groups for EMD, chi-square and Jensen–Shannon. |

## Returns

- Vector metrics: `{"sample_DR": DataFrame}`, a symmetric sample × sample matrix scaled to [0, 1]. Raises `ValueError` if `uns['X_DR_sample']` is missing.
- `"EMD"`: `{"EMD": DataFrame}`.
- `"chi_square"`, `"jensen_shannon"`, unknown methods: `None` (results are only written to disk).

## Output files

| Method | Folder | Files |
| --- | --- | --- |
| vector metric | `{output_dir}/{method}/` | `sample_DR_coordinates.csv`, `distance_matrix_sample_DR.csv`, `sample_distance_sample_DR_heatmap.pdf`, `distance_statistics_summary_{method}.csv`, `distance_check_results_sample_DR_{method}.txt` (group check) |
| `"EMD"` | `{output_dir}/EMD_distance/` | `distance_matrix_EMD.csv`, `cell_type_proportions.csv`, `cell_type_centroids.csv`, `distance_matrix_EMD_heatmap.pdf`, `distance_check_results_EMD_EMD.txt` |
| `"chi_square"` / `"jensen_shannon"` | `{output_dir}/chi_square/` or `{output_dir}/jensen_shannon/` | `distance_statistics_summary_<method>.csv`, and under `proportion_DR_distance/`: `proportion_DR_coordinates.csv`, `distance_matrix_proportion_DR.csv`, `cell_type_abundances.pdf`, `sample_distance_proportion_DR_heatmap.pdf`, `distance_check_results_proportion_DR_<method>.txt` |

!!! bug "`chi_square` and `jensen_shannon` fail on 0.5.0"
    In 0.5.0 both methods stop with `NameError: name 'plot_cell_type_abundances' is not defined`, right after writing `proportion_DR_coordinates.csv` and `distance_matrix_proportion_DR.csv`, and their group check is skipped because `grouping_columns` is not passed on. In the wrapper, listing either method in `*_sample_distance_methods` stops the modality's run. Both are fixed in 0.6.0 (on the demo, with `grouping_columns=["sev.level"]`, the group-check scores were 1.215 for chi-square and 1.122 for Jensen–Shannon). On 0.5.0, use a vector metric or `"EMD"`.

## Usage

Continuing from the RNA chain:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_distance.sample_distance import sample_distance

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

cos = sample_distance(pseudo, "api_demo/rna/Sample_distance", method="cosine",
                      data_type="RNA", grouping_columns=["sev.level"])
print(cos["sample_DR"].round(2))

emd = sample_distance(pseudo, "api_demo/rna/Sample_distance", method="EMD",
                      data_type="RNA", grouping_columns=["sev.level"],
                      cell_adata=adata, pseudobulk_adata=pseudo)
```

**Writes** → under `api_demo/rna/Sample_distance/`:

- `cosine/`: `distance_matrix_sample_DR.csv`, `sample_DR_coordinates.csv`, `sample_distance_sample_DR_heatmap.pdf`, `distance_statistics_summary_cosine.csv`, `distance_check_results_sample_DR_cosine.txt`
- `EMD_distance/`: `distance_matrix_EMD.csv`, `cell_type_proportions.csv`, `cell_type_centroids.csv`, `distance_matrix_EMD_heatmap.pdf`, `distance_check_results_EMD_EMD.txt`

The group check prints `Score for sample_DR (cosine): 1.385127` and `Score for EMD (EMD): 1.039446`: severe and healthy samples are further apart than samples within a group, clearly so in the embedding and barely in cell-type proportions alone. Both calls together take about 2 s (2.0 GB).
