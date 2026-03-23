# `preprocess_linux` (ATAC)

**Module**: `code/preparation/atac_preprocess_gpu.py`

```python
preprocess_linux(
    h5ad_path,
    sample_meta_path,
    output_dir,
    cell_meta_path=None,
    sample_column="sample",
    sample_level_batch_key="batch",
    cell_embedding_num_PCs=50,
    num_harmony_iterations=30,
    num_cell_hvfs=50000,
    min_cells=1,
    min_features=2000,
    max_features=15000,
    min_cells_per_sample=1,
    exclude_features=None,
    cell_level_batch_key=None,
    doublet_detection=True,
    tfidf_scale_factor=1e4,
    log_transform=True,
    drop_first_lsi=True,
    verbose=True,
)
```

End-to-end ATAC preprocessing on Linux with TF-IDF/LSI workflow.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `h5ad_path` | required | Input ATAC AnnData path. |
| `sample_meta_path` | required | Sample metadata CSV for merge into `obs`. |
| `output_dir` | required | Output root; writes into `preprocess/`. |
| `cell_meta_path` | `None` | Optional cell metadata CSV. |
| `sample_column` | `"sample"` | Sample ID column in `obs`. |
| `sample_level_batch_key` | `"batch"` | Sample-level batch column(s). |
| `cell_embedding_num_PCs` | `50` | LSI components used for cell embedding. |
| `num_harmony_iterations` | `30` | Harmony iteration cap. |
| `num_cell_hvfs` | `50000` | Number of highly variable features. |
| `min_cells` | `1` | Min cells with feature threshold. |
| `min_features` | `2000` | Minimum features per cell. |
| `max_features` | `15000` | Maximum features per cell. |
| `min_cells_per_sample` | `1` | Drop samples below this cell count. |
| `exclude_features` | `None` | Feature list to remove before modeling. |
| `cell_level_batch_key` | `None` | Cell-level Harmony batch keys. |
| `doublet_detection` | `True` | Enable scrublet-style doublet filtering. |
| `tfidf_scale_factor` | `1e4` | TF-IDF normalization scale factor. |
| `log_transform` | `True` | Apply log1p after TF-IDF. |
| `drop_first_lsi` | `True` | Drop first LSI component. |
| `verbose` | `True` | Print progress logs. |

## Returns

- `adata_cluster`: batch-corrected clustering object
- `adata_sample_diff`: sample-level analysis object

