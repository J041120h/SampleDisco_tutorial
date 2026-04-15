# `calculate_sample_embedding`

Aggregates single-cell data into a sample-level pseudobulk AnnData and computes two complementary sample embeddings:

1. **Expression embedding** (`X_DR_expression`) — PCA on the HVG-filtered pseudobulk matrix, optionally followed by Harmony and ComBat corrections. Captures the cell-type-weighted transcriptional signal of each sample.
2. **Proportion embedding** (`X_DR_proportion`) — PCA on the vector of cell-type proportions per sample, optionally Harmony-corrected. Captures how each sample differs in cellular composition.

Works for both RNA and ATAC: set `atac=True` to switch to TF-IDF/LSI-aware aggregation. GPU acceleration is enabled via `use_gpu=True`, which routes through `compute_pseudobulk_adata_linux`.

**Source:** `sample_embedding/calculate_sample_embedding.py:17`

## Signature

```python
def calculate_sample_embedding(
    adata: sc.AnnData,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    batch_col: Optional[Union[str, List[str]]] = None,
    output_dir: str = "./",
    sample_hvg_number: int = 2000,
    combat_timeout: Optional[float] = None,
    n_expression_components: int = 10,
    n_proportion_components: int = 10,
    harmony_for_proportion: bool = True,
    preserve_cols_in_sample_embedding: Optional[Union[str, List[str]]] = None,
    use_gpu: bool = False,
    atac: bool = False,
    save: bool = True,
    verbose: bool = True,
) -> Tuple[Dict, sc.AnnData]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level AnnData with a `cell_type` column in `.obs`. |
| `sample_col` | str | `"sample"` | Column in `.obs` identifying samples. |
| `celltype_col` | str | `"cell_type"` | Column in `.obs` identifying cell types. |
| `batch_col` | str or list, optional | `None` | Sample-level batch column(s); enables ComBat + Harmony on the expression embedding. |
| `output_dir` | str | `"./"` | Writes to `{output_dir}/pseudobulk/` and `{output_dir}/embeddings/`. |
| `sample_hvg_number` | int | `2000` | Number of HVGs per cell type used when building the pseudobulk matrix. |
| `combat_timeout` | float, optional | `None` | ComBat time limit in seconds. Defaults to 20s (GPU) or 1800s (CPU). |
| `n_expression_components` | int | `10` | PCA components for the expression embedding. |
| `n_proportion_components` | int | `10` | PCA components for the proportion embedding. |
| `harmony_for_proportion` | bool | `True` | Apply Harmony to the proportion embedding. |
| `preserve_cols_in_sample_embedding` | str or list, optional | `None` | Metadata columns to carry through to the pseudobulk `.obs`. |
| `use_gpu` | bool | `False` | Use the GPU implementation (`compute_pseudobulk_adata_linux`). |
| `atac` | bool | `False` | Switch to TF-IDF/LSI-aware aggregation for ATAC data. |
| `save` | bool | `True` | Write the pseudobulk AnnData and CSV exports to disk. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`Tuple[Dict, AnnData]` — `(pseudobulk_result_dict, sample_embedding_adata)`:

- `pseudobulk_result_dict` has keys `"cell_expression_corrected"` and `"cell_proportion"` (pandas DataFrames of the sample × feature matrices before DR).
- `sample_embedding_adata` is the sample-level AnnData with `X_DR_expression` and `X_DR_proportion` in `.obsm`.

## Output files

- `{output_dir}/pseudobulk/pseudobulk_sample.h5ad`
- `{output_dir}/embeddings/X_DR_expression.csv`
- `{output_dir}/embeddings/X_DR_proportion.csv`

## Usage

```python
from genodistance.sample_embedding import calculate_sample_embedding

pseudo_dict, pseudo_adata = calculate_sample_embedding(
    adata=adata_sample,
    sample_col="sample",
    celltype_col="cell_type",
    output_dir="/results/rna",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    use_gpu=True,
    atac=False,
)
```
