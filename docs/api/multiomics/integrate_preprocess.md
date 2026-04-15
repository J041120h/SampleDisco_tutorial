# `integrate_preprocess`

Post-integration QC for the GLUE-merged AnnData. Reads the integrated object produced by [`multiomics_preparation`](multiomics_preparation.md), merges per-modality sample metadata, applies filtering thresholds (minimum cells per sample, minimum genes per cell, mitochondrial cutoff, optional doublet detection), and writes the cleaned object back to disk. Duplicate sample IDs across modalities are deduplicated automatically so each sample shows up once per modality.

**Source:** `preparation/multi_omics_preprocess.py:147`

## Signature

```python
def integrate_preprocess(
    output_dir,
    h5ad_path=None,
    sample_column="sample",
    modality_col="modality",
    min_cells_sample=1,
    min_cell_gene=10,
    min_features=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    doublet=True,
    verbose=True,
    original_sample_col=None,
    rna_sample_meta_file=None,
    atac_sample_meta_file=None,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `output_dir` | str | — | Workspace directory; reads/writes `{output_dir}/preprocess/adata_sample.h5ad`. |
| `h5ad_path` | str, optional | `None` | Override the input h5ad. Defaults to `{output_dir}/preprocess/adata_sample.h5ad`. |
| `sample_column` | str | `"sample"` | Sample identifier column. |
| `modality_col` | str | `"modality"` | Column in `.obs` marking RNA vs ATAC. |
| `min_cells_sample` | int | `1` | Drop samples with fewer than this many cells. Bumped to 30 automatically if `doublet=True`. |
| `min_cell_gene` | int | `10` | Drop genes expressed in fewer than this many cells. |
| `min_features` | int | `500` | Drop cells with fewer than this many expressed features. |
| `pct_mito_cutoff` | float | `20` | Drop cells whose mitochondrial percentage exceeds this. |
| `exclude_genes` | list, optional | `None` | Genes/features to drop before filtering. |
| `doublet` | bool | `True` | Run doublet detection. |
| `verbose` | bool | `True` | Print progress. |
| `original_sample_col` | str, optional | `None` | Column storing the pre-deduplication sample IDs (defaults to `f"original_{sample_column}"`). |
| `rna_sample_meta_file`, `atac_sample_meta_file` | str, optional | `None` | Per-modality sample metadata CSVs merged in during this stage. |

## Returns

None. Updates `{output_dir}/preprocess/adata_sample.h5ad` in place.

## Output files

- `{output_dir}/preprocess/adata_sample.h5ad` (overwritten).

## Usage

```python
from genodistance.preparation import integrate_preprocess

integrate_preprocess(
    output_dir="/results/multiomics",
    sample_column="sample",
    modality_col="modality",
    min_cells_sample=1,
    min_features=500,
    pct_mito_cutoff=20,
    doublet=True,
)
```
