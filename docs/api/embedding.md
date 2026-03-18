# API Reference: Sample Embedding

These functions create the sample-level representations that power downstream distance, trajectory, and clustering analyses.

## `calculate_sample_embedding(...)`

```python
calculate_sample_embedding(
    adata: sc.AnnData,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    batch_col: Optional[Union[str, List[str]]] = None,
    output_dir: str = "./",
    sample_hvg_number: int = 2000,
    n_expression_components: int = 10,
    n_proportion_components: int = 10,
    harmony_for_proportion: bool = True,
    use_gpu: bool = False,
    atac: bool = False,
    save: bool = True,
    verbose: bool = True,
) -> Tuple[Dict, sc.AnnData]
```

Builds sample-level embeddings from single-cell data by computing pseudobulk summaries and then reducing expression-derived and proportion-derived representations.

**Key parameters**

- `adata`: input cell-level AnnData.
- `sample_col`, `celltype_col`: grouping columns needed for pseudobulk generation.
- `batch_col`: optional batch-correction columns.
- `sample_hvg_number`: number of highly variable features used per cell type.
- `n_expression_components`, `n_proportion_components`: final embedding dimensions.
- `atac=True`: switch to the ATAC-specific TF-IDF/LSI path.

**Returns**

- A pseudobulk result dictionary with expression and proportion components.
- A sample-level AnnData object with reduced coordinates stored in `.obsm`.

**Example**

```python
from sample_embedding.calculate_sample_embedding import calculate_sample_embedding

pseudobulk, sample_adata = calculate_sample_embedding(
    adata=adata_sample,
    sample_col="sample",
    celltype_col="cell_type",
    output_dir="results/rna",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
)
```

## `calculate_multiomics_sample_embedding(...)`

```python
calculate_multiomics_sample_embedding(
    adata: sc.AnnData,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    batch_col: Optional[Union[str, List[str]]] = None,
    output_dir: str = "./",
    sample_hvg_number: int = 2000,
    n_expression_components: int = 10,
    n_proportion_components: int = 10,
    harmony_for_proportion: bool = True,
    hvg_modality: str = "RNA",
    modality_col: str = "modality",
    ...
) -> Tuple[Dict, sc.AnnData]
```

Equivalent to `calculate_sample_embedding(...)`, but designed for integrated multi-omics objects with a modality label.

**Key parameters**

- `hvg_modality`: modality used for feature selection.
- `modality_col`: observation column identifying RNA vs ATAC cells.
- `preserve_cols_in_sample_embedding`: metadata columns that should survive correction steps.

**Returns**

- A pseudobulk result dictionary and a sample-level embedding AnnData.

## `find_optimal_cell_resolution(...)`

```python
find_optimal_cell_resolution(
    adata_cell: AnnData,
    adata_sample: AnnData,
    output_dir: str,
    column: str,
    modality: Literal["rna", "atac"] = "rna",
    coarse_start: float = 0.1,
    coarse_end: float = 1.0,
    coarse_step: float = 0.1,
    fine_range: float = 0.02,
    fine_step: float = 0.01,
    ...
) -> Tuple[float, pd.DataFrame]
```

Searches for an optimal cell-type clustering resolution by measuring concordance between cell-level structure and sample-level trajectory alignment.

**Returns**

- The selected resolution.
- A DataFrame summarizing resolution candidates and their scores.

## `replace_optimal_dimension_reduction(...)`

Multi-omics helper used during integrated resolution selection to swap in the preferred expression or proportion representation after optimization.

## See also

- [Preparation API](preparation.md)
- [Downstream Analysis API](downstream.md)
- [Tutorial 3: Multi-omics Integration](../tutorials/tutorial_multiomics.md)
