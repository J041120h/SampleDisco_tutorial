# `preprocess` (RNA)

End-to-end preprocessing for scRNA-seq on the CPU code path (a drop-in GPU variant lives in `rna_preprocess_gpu.preprocess_gpu`). The function ingests a cell-level `.h5ad`, merges sample and cell metadata, applies quality-control filters, normalizes and log-transforms, flags highly variable genes, runs PCA, and corrects batch effects with **two Harmony passes**. It returns a single AnnData carrying both cell embeddings: `Z_comp` (sample-removed, used for clustering and the composition blocks) and `Z_rmd` (sample-preserved, used for the RMD displacement block). Raw counts are preserved in `.layers['counts']` for downstream sample-level differential analysis, and the result is written atomically to disk so the pipeline can resume mid-way by reading the file back in.

**Source:** `preparation/rna_preprocess_cpu.py:167` (GPU: `preparation/rna_preprocess_gpu.py`)

## Signature

```python
def preprocess(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
    cell_embedding_num_PCs=20,
    num_harmony_iterations=30,
    num_cell_hvgs=2000,
    min_cells=500,
    min_genes=500,
    pct_mito_cutoff=20,
    exclude_genes=None,
    cell_level_batch_key=None,
    verbose=True,
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `h5ad_path` | str | — | Path to the cell-level `.h5ad` input. |
| `sample_meta_path` | str | — | CSV with one row per sample, keyed by `sample_column`. Pass `None` to skip the merge. |
| `output_dir` | str | — | Root output directory; the function writes to `{output_dir}/preprocess/`. |
| `sample_column` | str | `"sample"` | Column in `.obs` that identifies samples. Inferred from `obs_names` (split on `:`) if absent. |
| `cell_meta_path` | str, optional | `None` | Optional per-cell metadata CSV (indexed by `barcode`) merged into `adata.obs`. |
| `sample_level_batch_key` | str or list, optional | `None` | Sample-level batch column(s) validated against `.obs`. Set to `None` to skip. |
| `cell_embedding_num_PCs` | int | `20` | Number of PCs used in the cell-level embedding that Harmony corrects. |
| `num_harmony_iterations` | int | `30` | Max Harmony iterations. |
| `num_cell_hvgs` | int | `2000` | Number of highly variable genes (Seurat v3) flagged for the embedding. |
| `min_cells` | int | `500` | Drop genes expressed in fewer than this many cells; also the minimum cells-per-sample threshold. |
| `min_genes` | int | `500` | Drop cells with fewer than this many expressed genes. |
| `pct_mito_cutoff` | float | `20` | Drop cells whose mitochondrial percentage exceeds this. |
| `exclude_genes` | list, optional | `None` | Genes to drop before HVG selection (e.g. ribosomal). |
| `cell_level_batch_key` | list, optional | `None` | Cell-level batch keys for Harmony. The sample id is always added to the pass-1 (sample-removed) key automatically and always excluded from the pass-2 (sample-preserved) key. |
| `verbose` | bool | `True` | Print progress. |

## Returns

`AnnData` — a single preprocessed cell-level object with:

- `.X` — normalized + log1p expression (all genes retained).
- `.layers['counts']` — original raw counts (for pseudobulk).
- `.var['highly_variable']` — HVG flag (no subsetting).
- `.obsm['X_pca']` — PCA on the HVG subset.
- `.obsm['Z_comp']` — Harmony pass 1, **sample-removed** (composition view); input to [`cell_types`](cell_types_linux.md). Also written under the legacy alias `Z_clust` so the file stays readable by SampleDisco 0.2.0.
- `.obsm['X_pca_rmd']` — PCA on the sample-naive HVG subset; the basis Harmony pass 2 corrects.
- `.obsm['Z_rmd']` — Harmony pass 2, **sample-preserved** (RMD/displacement view).

Older h5ads carry earlier names for these two views (`Z_clust`/`Z_cmd`, or `X_pca_harmony`/`X_pca_harmony_nosamp`). They are still read automatically — see [Legacy embedding key names](../../legacy_keys.md).

!!! note "Single-key pipeline"
    This is the current single-output path. The old API returned a 2-tuple `(adata_cluster, adata_sample)` and wrote two files; that has been replaced by one AnnData and one file.

## Output files

- `{output_dir}/preprocess/adata_preprocessed.h5ad`

## Usage

```python
from sampledisco.preparation.rna_preprocess_cpu import preprocess

adata = preprocess(
    h5ad_path="data/test_RNA.h5ad",
    sample_meta_path=None,
    output_dir="sampledisco_demo_output/rna",
    sample_column="sample",
    cell_level_batch_key=None,
)
```

!!! tip "Prefer the config-driven wrapper"
    Most users should not call `preprocess` directly. The supported entry point runs preprocessing, cell-type clustering, sample embedding, and all downstream analysis from one YAML config:

    ```bash
    sampledisco -m complex --config config.yaml
    ```

    The GPU path (RAPIDS-accelerated normalization and Harmony) activates automatically when the RAPIDS stack is importable; otherwise SampleDisco falls back cleanly to this CPU implementation.
