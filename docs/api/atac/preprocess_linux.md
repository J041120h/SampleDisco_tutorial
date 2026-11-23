# `preprocess` (ATAC)

End-to-end preprocessing for scATAC-seq. It reads a cell × peak `.h5ad`, merges optional metadata, filters peaks, cells and samples, optionally removes Scrublet doublets, applies TF-IDF (+ `log1p`), and builds the two cell embeddings SampleDisco needs:

- `Z_comp`, **sample removed**: LSI on a sample-aware highly-variable-feature (HVF) set, then Harmony with the sample (plus `cell_level_batch_key`) as covariates.
- `Z_rmd`, **sample preserved**: LSI on a second, sample-naive HVF set, then Harmony on `cell_level_batch_key` only. With no `cell_level_batch_key`, this pass is skipped and `Z_rmd` is that LSI as is.

Raw counts are kept in `.layers['counts']`. The result is written to `{output_dir}/preprocess/adata_preprocessed.h5ad` and returned.

**Import:** `from sampledisco.preparation.atac_preprocess_cpu import preprocess`

**Source:** `preparation/atac_preprocess_cpu.py:184`. GPU variant: `preprocess_gpu` in `preparation/atac_preprocess_gpu.py:167`.

## Signature

```python
def preprocess(
    h5ad_path,
    sample_meta_path,
    output_dir,
    sample_column="sample",
    cell_meta_path=None,
    sample_level_batch_key=None,
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
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `h5ad_path` | str | — | Cell × peak `.h5ad` with raw counts. |
| `sample_meta_path` | str or None | — (required positional) | Per-sample metadata CSV merged into `.obs`; `None` to skip. |
| `output_dir` | str | — | The file goes to `{output_dir}/preprocess/`. |
| `sample_column` | str | `"sample"` | Sample column; inferred from `obs_names` (text before the first `:`) if missing. |
| `cell_meta_path` | str or None | `None` | Per-cell metadata CSV indexed by a `barcode` column. |
| `sample_level_batch_key` | str, list or None | `None` | Only checked to exist in `.obs`; it does not reach Harmony. Sample-level batch correction happens in [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) via `batch_col`. |
| `cell_embedding_num_PCs` | int | `50` | LSI components computed. With `drop_first_lsi=True` the embeddings keep `cell_embedding_num_PCs − 1` (49 by default). |
| `num_harmony_iterations` | int | `30` | Maximum Harmony iterations. |
| `num_cell_hvfs` | int | `50000` | Features per HVF set (Seurat v3 on the TF-IDF matrix); flagged, not subset. |
| `min_cells` | int | `1` | Drop peaks detected in fewer cells. |
| `min_features` / `max_features` | int | `2000` / `15000` | Keep cells whose detected-peak count is within these bounds. |
| `min_cells_per_sample` | int | `1` | Drop samples with fewer surviving cells. |
| `exclude_features` | list or None | `None` | Peak IDs removed after doublet removal. |
| `cell_level_batch_key` | list or None | `None` | Cell-level batch columns. Pass 1 (`Z_comp`) uses these plus `sample_column`; pass 2 (`Z_rmd`) uses these without `sample_column`. |
| `doublet_detection` | bool | `True` | Run Scrublet per sample and drop predicted doublets (skipped with a message if Scrublet fails). |
| `tfidf_scale_factor` | float | `1e4` | Scale factor in TF-IDF. |
| `log_transform` | bool | `True` | `log1p` after TF-IDF. |
| `drop_first_lsi` | bool | `True` | Drop LSI component 1, which mostly tracks sequencing depth. |
| `verbose` | bool | `True` | Print progress. |

After these filters, peaks detected in fewer than 0.1% of the remaining cells are dropped.

## Returns

`AnnData`, also written to disk. `.X` holds TF-IDF (+ `log1p`) values, `.layers['counts']` the raw counts, and `.var['highly_variable']` / `.var['highly_variable_rmd']` the two HVF sets. `.obsm` carries:

- `X_lsi`: LSI on the sample-aware HVFs (first component dropped by default). Its presence is how [`cell_types`](../rna/cell_types_linux.md) recognizes ATAC.
- `Z_comp`: Harmony pass 1, sample removed. Also written as the 0.2.0 alias `Z_clust`.
- `X_lsi_rmd`: LSI on the sample-naive HVFs.
- `Z_rmd`: Harmony pass 2, sample preserved (or `X_lsi_rmd` itself without `cell_level_batch_key`).

Older files use earlier names for these views; see [Legacy embedding key names](../../legacy_keys.md).

## Output files

- `{output_dir}/preprocess/adata_preprocessed.h5ad`

## Usage

The ATAC demo has no batch column, so no batch key is passed:

```python
from sampledisco.preparation.atac_preprocess_cpu import preprocess

adata_atac = preprocess(
    h5ad_path="data/test_ATAC.h5ad",
    sample_meta_path=None,
    output_dir="api_demo/atac",
    sample_column="sample",
)
print(adata_atac)
```

**Writes** → `api_demo/atac/preprocess/adata_preprocessed.h5ad`. From the log: QC keeps 25,196 of 29,021 cells, Scrublet removes 3 doublets, and 228,586 peaks remain; pass 1 uses `sample` as the only covariate and pass 2 is skipped (`no extra batch covariate → using raw X_lsi_rmd`). The returned object:

```text
AnnData object with n_obs × n_vars = 25193 × 228586
    ...
    obsm: 'X_lsi', 'Z_comp', 'Z_clust', 'X_lsi_rmd', 'Z_rmd'
    varm: 'LSI'
    layers: 'counts'
```

`Z_comp`, `X_lsi_rmd` and `Z_rmd` have 49 columns (50 LSI components minus the first). 5.5–6.5 minutes and 8.3 GB peak RAM on 8 cores (two runs).

Next step: [`cell_types`](../rna/cell_types_linux.md) (the ATAC example there), then [`compute_sample_embedding`](../shared/calculate_sample_embedding.md).

!!! note "GPU"
    `preprocess_gpu` runs only when called directly or when the wrapper config sets `use_gpu: true` (with a CPU fallback if RAPIDS does not import). It takes the same arguments and defaults, but `cell_meta_path` comes before `sample_column` positionally. See [GPU acceleration](../../tutorials/gpu.md).

!!! note "Prefer the config-driven wrapper"
    `sampledisco -m complex --config <config.yaml>` with `run_atac_pipeline: true` runs ATAC preprocessing → [`cell_types`](../rna/cell_types_linux.md) → sample embedding → downstream analysis.
