# `preprocess_linux` (RNA)

**Module**: `code/preparation/rna_preprocess_gpu.py`

```python
preprocess_linux(
    h5ad_path,
    sample_meta_path,
    output_dir,
    cell_meta_path=None,
    sample_column="sample",
    sample_level_batch_key="batch",
    cell_embedding_num_PCs=20,
    num_harmony_iterations=30,
    num_cell_hvgs=2000,
    min_cells=500,
    min_genes=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    cell_level_batch_key=None,
    verbose=True,
)
```

End-to-end RNA preprocessing on Linux with GPU-accelerated downstream steps.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `h5ad_path` | required | Input RNA AnnData file path. |
| `sample_meta_path` | required | Sample metadata CSV path for merging into `obs`. |
| `output_dir` | required | Output root; function writes into `preprocess/`. |
| `cell_meta_path` | `None` | Optional cell metadata CSV for `obs` augmentation. |
| `sample_column` | `"sample"` | Sample ID column name in `obs`. |
| `sample_level_batch_key` | `"batch"` | Sample-level batch column(s). |
| `cell_embedding_num_PCs` | `20` | PCA components for cell-level embedding. |
| `num_harmony_iterations` | `30` | Harmony iteration limit. |
| `num_cell_hvgs` | `2000` | Number of highly variable genes. |
| `min_cells` | `500` | Min cells per gene filter threshold. |
| `min_genes` | `500` | Min genes per cell filter threshold. |
| `pct_mito_cutoff` | `20` | Max mitochondrial percentage cutoff. |
| `exclude_genes` | `None` | Gene list to exclude before processing. |
| `cell_level_batch_key` | `None` | Cell-level batch keys for Harmony. |
| `verbose` | `True` | Print progress messages. |

## Returns

- `adata_cluster`: batch-corrected object for clustering
- `adata_sample_diff`: minimally processed object for sample-level analysis

