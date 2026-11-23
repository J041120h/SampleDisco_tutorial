# `build_sample_adata`

Turns the sample embedding stored on a cell-level AnnData into a small **sample-level AnnData** with one row per unit (sample, or `<sample>_<modality>` in multi-omics). This is the `pseudo_adata` / `pseudobulk_adata` that the [downstream](../downstream.md) functions take. The wrapper builds it the same way before calling [`downstream_analysis`](../pipeline/downstream_analysis.md).

**Import:** `from sampledisco.sample_embedding.sample_embedding import build_sample_adata`

**Source:** `sample_embedding/sample_embedding.py:64`

## Signature

```python
def build_sample_adata(
    adata,
    sample_col: str = "sample",
    modality_col: Optional[str] = None,
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level object with `uns['X_DR_sample']`, from [`compute_sample_embedding`](calculate_sample_embedding.md) or [`run_autotune`](run_autotune.md), or reloaded from an h5ad that stored it. A `KeyError` is raised if the key is missing. |
| `sample_col` | str | `"sample"` | Sample column in `adata.obs`. |
| `modality_col` | str or None | `None` | **Set it for multi-omics** (e.g. `"modality"`), so cells are grouped into the same `<sample>_<modality>` units as the embedding. Without it, the multi-omics unit ids do not match and the metadata comes out as `NaN`. |

## Returns

A new `AnnData` with one row per embedding unit, in the embedding's order:

- `.X`: the embedding values (float32, units × PCs).
- `.obsm['X_DR_sample']`: the same matrix; `.uns['X_DR_sample']`: a copy of the `DataFrame`.
- `.uns['sample_embedding_params']`: copied when present.
- `.obs`: indexed by unit id (index name `sample`), with the `adata.obs` columns that are **constant within every unit** (missing values ignored), holding that value. Columns that vary within a unit (per-cell values such as `n_genes`, `total_counts` or the cell types) are dropped and listed in one printed line, `[build_sample_adata] dropped N obs column(s) that vary within a sample (per-cell values): [...]`.

## Output files

None; the object lives in memory.

## Usage

Continuing from the [`compute_sample_embedding`](calculate_sample_embedding.md) example. This is also the setup block that every RNA downstream example starts with:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")  # cells
pseudo = build_sample_adata(adata, sample_col="sample")                # one row per sample
print(pseudo)
print(pseudo.obs[["sev.level", "batch"]])
```

Output (under a second):

```text
[build_sample_adata] dropped 10 obs column(s) that vary within a sample (per-cell values): ['celltype', 'old', 'large', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'cell_type_original']
AnnData object with n_obs × n_vars = 8 × 7
    obs: 'sample', 'batch', 'sev.level'
    uns: 'X_DR_sample', 'sample_embedding_params'
    obsm: 'X_DR_sample'
             sev.level  batch
sample
CoV-3-Wilk           4   Wilk
CoV-34-Aruna         4  Aruna
CoV-36-Aruna         4  Aruna
CoV-4-Wilk           4   Wilk
HD-1-Wilk            1   Wilk
HD-20-Aruna          1  Aruna
HD-21-Aruna          1  Aruna
HD-4-Wilk            1   Wilk
```

For multi-omics, see the example on [`compute_sample_embedding` (multi-omics)](../multiomics/calculate_multiomics_sample_embedding.md), which passes `modality_col="modality"`.
