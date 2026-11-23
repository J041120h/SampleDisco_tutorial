# `harmonize_xglue`

Derives the **sample-removed** view `Z_comp` from a sample-preserved joint embedding by one Harmony pass with the sample as covariate. scGLUE's `obsm['X_glue']` keeps per-sample variance, so it serves as `Z_rmd`; this function produces the `Z_comp` that cell typing and the composition blocks need. The multi-omics wrapper calls it after `multiomics_preparation` whenever the union has no composition key yet (i.e. unless `multiomics_run_glue_twice_for_sample_removal: true` already produced `Z_comp` from a second scGLUE run).

What it does:

1. Batch covariates: `[batch_col, sample_col]` if `batch_col` is in `.obs` with at least 2 levels, else `[sample_col]`.
2. Harmony ([`harmonize_embedding`](../utils/harmonize_embedding.md)) on `obsm[in_key]` → `obsm[out_key]` (float32). When `out_key` is `"Z_comp"`, the 0.2.0 alias `obsm['Z_clust']` is written too (overwriting an existing one).
3. If `obsm['Z_rmd']` does not exist yet, `obsm[in_key]` is copied to `obsm['Z_rmd']`.

It is a no-op (with a log line) when `in_key` is missing from `.obsm` or `sample_col` has fewer than 2 levels. Nothing is written to disk.

It works for any joint embedding, not only scGLUE: the [MaxFuse tutorial](../../tutorials/multiomics_maxfuse.md) uses it with `in_key="Z_rmd"` on a MaxFuse embedding.

**Import:** `from sampledisco.preparation.multi_omics_batch_correction import harmonize_xglue`

**Source:** `preparation/multi_omics_batch_correction.py:58`

## Signature

```python
def harmonize_xglue(
    adata: AnnData,
    *,
    sample_col: str,
    batch_col: Optional[str] = None,
    in_key: str = "X_glue",
    out_key: str = "Z_comp",
    use_gpu: bool = False,
    max_iter: int = 50,
    random_state: int = 0,
    verbose: bool = True,
) -> AnnData
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level object (e.g. the union from [`multiomics_preparation`](multiomics_preparation.md)). Modified in place. |
| `sample_col` | str | — (keyword, required) | Sample column; always a Harmony covariate. |
| `batch_col` | str or None | `None` | Optional batch column, added as a covariate when it has at least 2 levels. Wrapper key: `multiomics_batch_col`. |
| `in_key` | str | `"X_glue"` | Input `.obsm` key (the sample-preserved embedding). |
| `out_key` | str | `"Z_comp"` | Output `.obsm` key. |
| `use_gpu` | bool | `False` | Passed to harmony-pytorch (ignored by the harmonypy fallback). |
| `max_iter` | int | `50` | Maximum Harmony iterations. Wrapper key: `multiomics_harmonize_xglue_max_iter`. |
| `random_state` | int | `0` | Harmony seed. |
| `verbose` | bool | `True` | Print progress. |

## Returns

The same `AnnData` with `obsm[out_key]`, possibly `obsm['Z_clust']` and `obsm['Z_rmd']` added (see above).

## Usage

First step of the multi-omics chain. The demo file stores scGLUE's output as `X_glue`; recomputing `Z_comp` from it is what the wrapper does when no composition key exists:

```python
import os
import anndata as ad
from sampledisco.preparation.multi_omics_batch_correction import harmonize_xglue

mo = ad.read_h5ad("data/test_multiomics_integrated.h5ad")
mo = harmonize_xglue(mo, sample_col="sample", in_key="X_glue", out_key="Z_comp")
print(sorted(mo.obsm.keys()))

os.makedirs("api_demo/multiomics/preprocess", exist_ok=True)
mo.write_h5ad("api_demo/multiomics/preprocess/adata_sample.h5ad")
```

**Writes** → `api_demo/multiomics/preprocess/adata_sample.h5ad` (by the last line). Output:

```text
[xglue-harmony] Z_comp: 59,010 × 50 dims, batch_keys=['sample'], gpu=False, max_iter=50
[xglue-harmony] Z_comp: done in 10.3s (also aliased obsm['X_glue'] → obsm['Z_rmd'])
['X_glue', 'X_umap', 'Z_clust', 'Z_cmd', 'Z_comp', 'Z_rmd']
```

The demo file has no batch column, so only `sample` is a covariate. About 12 s and 1.1 GB.
