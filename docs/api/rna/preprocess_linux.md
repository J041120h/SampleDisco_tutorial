# `preprocess` (RNA)

End-to-end preprocessing for scRNA-seq. It reads a cell-level `.h5ad` with raw counts, merges optional cell and sample metadata, filters cells, genes and samples, normalizes and log-transforms, flags highly variable genes, runs PCA, and builds the two cell embeddings SampleDisco needs:

- `Z_comp`, **sample removed**: Harmony on PCA with the sample (plus `cell_level_batch_key`) as covariates. Used for cell typing and the composition blocks.
- `Z_rmd`, **sample preserved**: PCA on a second, sample-naive HVG set, then Harmony on `cell_level_batch_key` only. With no `cell_level_batch_key`, this second Harmony pass is skipped and `Z_rmd` is that PCA as is. Used for the RMD displacement block.

Raw counts are kept in `.layers['counts']`. The result is written to `{output_dir}/preprocess/adata_preprocessed.h5ad` (a plain `sc.write` of a cleaned copy, so a crash mid-write can leave a partial file) and returned.

**Import:** `from sampledisco.preparation.rna_preprocess_cpu import preprocess`

**Source:** `preparation/rna_preprocess_cpu.py:175`. GPU variant: `preprocess_gpu` in `preparation/rna_preprocess_gpu.py:168`.

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
| `h5ad_path` | str | — | Cell-level `.h5ad` with raw counts in `.X`. |
| `sample_meta_path` | str or None | — (required positional) | CSV with one row per sample, merged into `.obs` on `sample_column`. Pass `None` to skip. |
| `output_dir` | str | — | Root output directory; the file goes to `{output_dir}/preprocess/`. |
| `sample_column` | str | `"sample"` | Sample column in `.obs`. If missing, it is inferred from `obs_names` (the part before the first `:`). |
| `cell_meta_path` | str or None | `None` | Per-cell metadata CSV indexed by a `barcode` column, left-joined into `.obs`. |
| `sample_level_batch_key` | str, list or None | `None` | Only checked to exist in `.obs` (a `KeyError` lists missing columns). It does not enter Harmony here; sample-level batch correction happens in [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) via `batch_col`. |
| `cell_embedding_num_PCs` | int | `20` | PCA components for both embeddings. |
| `num_harmony_iterations` | int | `30` | Maximum Harmony iterations. |
| `num_cell_hvgs` | int | `2000` | HVGs (Seurat v3) per HVG set. Genes are flagged in `.var['highly_variable']`, not subset. |
| `min_cells` | int | `500` | Drop genes detected in fewer cells, and drop samples with fewer cells. |
| `min_genes` | int | `500` | Drop cells with fewer detected genes. |
| `pct_mito_cutoff` | float | `20` | Drop cells with `pct_counts_mt` at or above this (genes starting `MT-`/`mt-`). |
| `exclude_genes` | list or None | `None` | Genes removed before HVG selection. `MT-` genes are always removed after the mito filter. |
| `cell_level_batch_key` | list or None | `None` | Cell-level batch columns for Harmony. Pass 1 (`Z_comp`) uses these **plus** `sample_column`; pass 2 (`Z_rmd`) uses these **without** `sample_column`. |
| `verbose` | bool | `True` | Print progress. |

After the filters above, genes detected in fewer than 0.1% of the remaining cells are also dropped.

## Returns

`AnnData`, also written to disk, with:

- `.X`: normalized (10,000 counts per cell) + `log1p` expression, all retained genes.
- `.layers['counts']`: raw counts.
- `.var['highly_variable']` (= `highly_variable_clust`, the sample-aware HVG set) and `.var['highly_variable_rmd']` (the sample-naive set).
- `.obsm['X_pca']`: PCA on the sample-aware HVGs.
- `.obsm['Z_comp']`: Harmony pass 1, sample removed. Also written as the 0.2.0 alias `Z_clust` (to be removed in 1.0).
- `.obsm['X_pca_rmd']`: PCA on the sample-naive HVGs.
- `.obsm['Z_rmd']`: Harmony pass 2 on `X_pca_rmd`, sample preserved (or `X_pca_rmd` itself without `cell_level_batch_key`).

Older files use earlier names for the two views; they are still read automatically. See [Legacy embedding key names](../utils/embedding_keys.md) and [`migrate_obsm_keys`](../utils/embedding_keys.md).

## Output files

- `{output_dir}/preprocess/adata_preprocessed.h5ad`

## Usage

```python
from sampledisco.preparation.rna_preprocess_cpu import preprocess

adata = preprocess(
    h5ad_path="data/test_RNA.h5ad",
    sample_meta_path=None,            # metadata is already in .obs
    output_dir="api_demo/rna",
    sample_column="sample",
    cell_level_batch_key=["batch"],   # Harmony on batch in both passes
)
print(adata)
```

**Writes** → `api_demo/rna/preprocess/adata_preprocessed.h5ad`. On the demo, QC keeps 29,680 of 29,989 cells and 9,936 of 14,881 genes (all 8 samples); the log shows the two Harmony passes (`batch, sample`, then `batch`) and the returned object:

```text
AnnData object with n_obs × n_vars = 29680 × 9936
    obs: 'celltype', 'sample', 'old', 'large', 'batch', 'sev.level', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'
    var: 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'highly_variable_clust', 'highly_variable_rmd'
    uns: 'X_name', 'hvg', 'log1p', 'pca'
    obsm: 'X_pca', 'Z_comp', 'Z_clust', 'X_pca_rmd', 'Z_rmd'
    varm: 'PCs'
    layers: 'counts'
```

45–55 s and 2.9 GB peak RAM on 8 CPU cores (two runs).

!!! note "GPU"
    The GPU variant runs only when you ask for it: call `preprocess_gpu` yourself, or set `use_gpu: true` in the wrapper config (the wrapper falls back to this CPU function if the RAPIDS modules do not import). Having RAPIDS installed is not enough. `preprocess_gpu` takes the same arguments and defaults, but `cell_meta_path` comes before `sample_column` positionally. See [GPU acceleration](../../tutorials/gpu.md).

!!! tip "Prefer the config-driven wrapper"
    `sampledisco -m complex --config config.yaml` with `run_rna_pipeline: true` runs preprocessing, cell typing, sample embedding and the downstream steps. The keys `rna_cell_level_batch_key`, `rna_sample_level_batch_col`, `rna_min_cells`, … map one-to-one onto the arguments above.
