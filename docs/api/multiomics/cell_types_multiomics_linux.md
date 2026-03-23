# `cell_types_multiomics_linux`

**Module**: `code/preparation/multi_omics_cell_type_gpu.py`

```python
cell_types_multiomics_linux(
    adata,
    modality_column="modality",
    rna_modality_value="RNA",
    atac_modality_value="ATAC",
    cell_type_column="cell_type",
    cluster_resolution=0.8,
    use_rep="X_glue",
    num_PCs=50,
    k_neighbors=15,
    transfer_metric="cosine",
    compute_umap=True,
    save=False,
    output_dir=None,
    defined_output_path=None,
    verbose=True,
    generate_plots=True,
)
```

Clusters RNA cells, then transfers labels to ATAC cells with GPU KNN/SNN strategy.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | Integrated AnnData containing both modalities. |
| `modality_column` | `"modality"` | Column identifying modality per cell. |
| `rna_modality_value` | `"RNA"` | Label value treated as RNA subset. |
| `atac_modality_value` | `"ATAC"` | Label value treated as ATAC subset. |
| `cell_type_column` | `"cell_type"` | Output label column name. |
| `cluster_resolution` | `0.8` | RNA Leiden resolution. |
| `use_rep` | `"X_glue"` | Embedding used for clustering and transfer. |
| `num_PCs` | `50` | Number of embedding components used. |
| `k_neighbors` | `15` | KNN neighborhood size for label transfer. |
| `transfer_metric` | `"cosine"` | Distance metric for neighbor search. |
| `compute_umap` | `True` | Compute UMAP embedding. |
| `save` | `False` | Save updated AnnData. |
| `output_dir` | `None` | Output directory when saving/plotting. |
| `defined_output_path` | `None` | Explicit file path for updated AnnData. |
| `verbose` | `True` | Print progress details. |
| `generate_plots` | `True` | Generate summary plots. |

## Returns

- Updated integrated `AnnData` with transferred cell type labels.

