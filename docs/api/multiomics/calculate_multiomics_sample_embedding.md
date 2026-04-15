# `calculate_multiomics_sample_embedding`

Multi-omics counterpart of [`calculate_sample_embedding`](../shared/calculate_sample_embedding.md). Aggregates the integrated RNA + ATAC AnnData into a sample-level pseudobulk and computes the two sample embeddings (`X_DR_expression`, `X_DR_proportion`). Uses a dedicated multi-omics pseudobulk routine that honors the `modality_col` label and selects HVGs from the modality specified by `hvg_modality` (typically RNA, since gene-level features are more interpretable).

**Source:** `sample_embedding/calculate_multiomics_sample_embedding.py:18`

## Signature

```python
def calculate_multiomics_sample_embedding(
    adata: sc.AnnData,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    batch_col: Optional[Union[str, List[str]]] = None,
    output_dir: str = "./",
    sample_hvg_number: int = 2000,
    n_expression_components: int = 10,
    n_proportion_components: int = 10,
    harmony_for_proportion: bool = True,
    preserve_cols_in_sample_embedding: Optional[Union[str, List[str]]] = None,
    use_gpu: bool = False,
    atac: bool = False,
    save: bool = True,
    verbose: bool = True,
    hvg_modality: str = "RNA",
    modality_col: str = "modality",
) -> Tuple[Dict, sc.AnnData]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Integrated cell-level AnnData with `cell_type` labels and a `modality_col`. |
| `sample_col` | str | `"sample"` | Sample identifier column. |
| `celltype_col` | str | `"cell_type"` | Cell-type column. |
| `batch_col` | str or list, optional | `None` | Sample-level batch. Modality tag is already added automatically. |
| `output_dir` | str | `"./"` | Writes `{output_dir}/pseudobulk/pseudobulk_sample.h5ad`. |
| `sample_hvg_number` | int | `2000` | HVGs per cell type for the pseudobulk. |
| `n_expression_components` | int | `10` | Components for `X_DR_expression`. |
| `n_proportion_components` | int | `10` | Components for `X_DR_proportion`. |
| `harmony_for_proportion` | bool | `True` | Harmony on the proportion embedding. |
| `preserve_cols_in_sample_embedding` | str or list, optional | `None` | `.obs` columns to carry into pseudobulk `.obs`. |
| `use_gpu` | bool | `False` | Use the GPU-accelerated pseudobulk. |
| `atac` | bool | `False` | LSI-style aggregation (rare for multi-omics; usually left `False`). |
| `save` | bool | `True` | Write the pseudobulk to disk. |
| `verbose` | bool | `True` | Print progress. |
| `hvg_modality` | str | `"RNA"` | Modality whose HVGs are used during pseudobulk feature selection. |
| `modality_col` | str | `"modality"` | Column identifying modality per cell. |

## Returns

`Tuple[Dict, AnnData]` — same shape as the single-modality version:

- `pseudobulk_result_dict["cell_expression_corrected"]`, `["cell_proportion"]`.
- `sample_embedding_adata` with `X_DR_expression` and `X_DR_proportion` in `.obsm`.

## Output files

- `{output_dir}/pseudobulk/pseudobulk_sample.h5ad`
- `{output_dir}/embeddings/*.csv`

## Usage

```python
from genodistance.sample_embedding import calculate_multiomics_sample_embedding

pseudo_dict, pseudo_adata = calculate_multiomics_sample_embedding(
    adata=adata_integrated,
    sample_col="sample",
    celltype_col="cell_type",
    modality_col="modality",
    hvg_modality="RNA",
    output_dir="/results/multiomics",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    use_gpu=True,
)
```
