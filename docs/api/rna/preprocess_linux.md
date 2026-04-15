# `preprocess_linux` (RNA)

End-to-end preprocessing for scRNA-seq on the GPU/Linux code path. The function ingests a cell-level `.h5ad`, merges sample and cell metadata, applies quality-control filters, normalizes and log-transforms, selects highly variable genes, runs PCA, and corrects batch effects with Harmony. It returns two AnnData objects: one Harmony-integrated copy destined for clustering, and one minimally processed copy that keeps raw counts for downstream sample-level differential analysis. Intermediate and final outputs are written atomically to disk so the pipeline can resume mid-way by reading the files back in.

**Source:** `preparation/rna_preprocess_gpu.py:108`

## Signature

```python
def preprocess_linux(
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
) -> Tuple[AnnData, AnnData]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `h5ad_path` | str | — | Path to the cell-level `.h5ad` input. |
| `sample_meta_path` | str | — | CSV with one row per sample, keyed by `sample_column`. |
| `output_dir` | str | — | Root output directory; the function writes to `{output_dir}/preprocess/`. |
| `cell_meta_path` | str, optional | `None` | Optional per-cell metadata CSV merged into `adata.obs`. |
| `sample_column` | str | `"sample"` | Column in `.obs` that identifies samples. |
| `sample_level_batch_key` | str or list, optional | `"batch"` | Sample-level batch column(s) used for sample-level Harmony. Set to `None` to skip. |
| `cell_embedding_num_PCs` | int | `20` | Number of PCs used in the cell-level embedding that Harmony corrects. |
| `num_harmony_iterations` | int | `30` | Max Harmony iterations. |
| `num_cell_hvgs` | int | `2000` | Number of highly variable genes for the clustering track. |
| `min_cells` | int | `500` | Drop genes expressed in fewer than this many cells. |
| `min_genes` | int | `500` | Drop cells with fewer than this many expressed genes. |
| `pct_mito_cutoff` | float | `20` | Drop cells whose mitochondrial percentage exceeds this. |
| `exclude_genes` | list, optional | `None` | Genes to drop before HVG selection (e.g. ribosomal). |
| `cell_level_batch_key` | list, optional | `None` | Cell-level batch keys for Harmony. Sample id is always included automatically. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`Tuple[AnnData, AnnData]` — `(adata_cluster, adata_sample)`:

- `adata_cluster` carries `X_pca_harmony` in `.obsm` and is the input to [`cell_types_linux`](cell_types_linux.md).
- `adata_sample` preserves raw counts (for pseudobulk) with `X_pca` in `.obsm`.

## Output files

- `{output_dir}/preprocess/adata_cell.h5ad`
- `{output_dir}/preprocess/adata_sample.h5ad`

## Usage

```python
from genodistance.preparation import preprocess_linux

adata_cluster, adata_sample = preprocess_linux(
    h5ad_path="/data/test_RNA.h5ad",
    sample_meta_path="/data/sample_meta.csv",
    output_dir="/results/rna",
    sample_column="sample",
    cell_level_batch_key=None,
)
```
