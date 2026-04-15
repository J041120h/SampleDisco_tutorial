# `cell_types_multiomics_linux`

Assigns cell-type labels on the integrated RNA + ATAC object. RNA cells are clustered with Leiden on the joint embedding (`X_glue` by default); ATAC cells then receive labels via cosine-metric k-nearest neighbors in the same embedding. The k-NN graph is Jaccard-weighted (SNN) and imbalanced clusters can be dropped automatically to keep label transfer robust. Optionally computes a UMAP on the joint space.

**Source:** `preparation/multi_omics_cell_type_gpu.py:19`

## Signature

```python
def cell_types_multiomics_linux(
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

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Integrated AnnData with the joint embedding in `.obsm[use_rep]`. |
| `modality_column` | str | `"modality"` | Column marking each cell's modality. |
| `rna_modality_value` | str | `"RNA"` | Value identifying RNA cells. |
| `atac_modality_value` | str | `"ATAC"` | Value identifying ATAC cells. |
| `cell_type_column` | str | `"cell_type"` | Destination column for labels. |
| `cluster_resolution` | float | `0.8` | Leiden resolution for RNA clustering. |
| `use_rep` | str | `"X_glue"` | Key in `.obsm` used for clustering and k-NN transfer. |
| `num_PCs` | int, optional | `50` | Trim the joint embedding to this many components before processing. |
| `k_neighbors` | int | `15` | Nearest RNA neighbors per ATAC cell during label transfer. |
| `transfer_metric` | str | `"cosine"` | Distance metric used in the k-NN search. |
| `compute_umap` | bool | `True` | Run UMAP on the joint embedding after labeling. |
| `save` | bool | `False` | Persist the updated AnnData to disk. |
| `output_dir` | str, optional | `None` | Directory for UMAP PNGs and the saved h5ad. |
| `defined_output_path` | str, optional | `None` | Override the exact h5ad save path. |
| `verbose` | bool | `True` | Print progress. |
| `generate_plots` | bool | `True` | Emit a modality-balance heatmap alongside the UMAP. |

## Returns

The input AnnData, now carrying `.obs[cell_type_column]` and (when requested) `X_umap` in `.obsm`.

## Output files

- `{output_dir}/preprocess/adata_sample.h5ad` (when `save=True`).
- UMAP and diagnostic PNGs under `{output_dir}/preprocess/`.

## Usage

```python
from genodistance.preparation import cell_types_multiomics_linux

adata = cell_types_multiomics_linux(
    adata=adata_integrated,
    modality_column="modality",
    cluster_resolution=0.8,
    use_rep="X_glue",
    num_PCs=50,
    k_neighbors=1,
    transfer_metric="cosine",
    compute_umap=True,
    save=True,
    output_dir="/results/multiomics",
)
```
