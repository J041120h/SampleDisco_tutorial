# `integrate_preprocess`

**Module**: `code/preparation/multi_omics_preprocess.py`

```python
integrate_preprocess(
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

Preprocess integrated RNA+ATAC AnnData before sample-level embedding.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `output_dir` | required | Output root for preprocess artifacts. |
| `h5ad_path` | `None` | Integrated input h5ad path (defaults to `{output_dir}/preprocess/adata_sample.h5ad`). |
| `sample_column` | `"sample"` | Sample ID column in `obs`. |
| `modality_col` | `"modality"` | Modality label column (`RNA`/`ATAC`). |
| `min_cells_sample` | `1` | Minimum cells per sample filter threshold. |
| `min_cell_gene` | `10` | Minimum genes detected per cell. |
| `min_features` | `500` | Minimum features detected per cell. |
| `pct_mito_cutoff` | `20` | Mitochondrial percentage filter threshold. |
| `exclude_genes` | `None` | Gene list to remove. |
| `doublet` | `True` | Enable doublet filtering logic. |
| `verbose` | `True` | Print progress logs. |
| `original_sample_col` | `None` | Backup column name to store original sample IDs. |
| `rna_sample_meta_file` | `None` | Optional RNA metadata merge file. |
| `atac_sample_meta_file` | `None` | Optional ATAC metadata merge file. |

## Returns

- Preprocessed integrated `AnnData`

