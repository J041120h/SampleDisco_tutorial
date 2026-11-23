# `save_embedding_to_h5ad`

Stores the sample embedding of an in-memory AnnData (`uns['X_DR_sample']` and `uns['sample_embedding_params']`) in an `.h5ad` file, normally the one the object was read from.

- If the file holds the same object as memory (identical `obs` and `var`; the same `obsm`/`layers`/`obsp`/`varm` keys; `X` and every matrix with the same shape, dtype and a fixed sample of values; the same `uns` key names), only the two `uns` entries are replaced, in place. Each entry is written under a temporary key and then moved, so an interrupted write leaves the old or the new value.
- Otherwise the whole file is rewritten from memory with `sc.write`.

[`compute_sample_embedding`](calculate_sample_embedding.md) and [`run_autotune`](run_autotune.md) call it themselves, but only for `{output_dir}/preprocess/adata_preprocessed.h5ad`. Call it yourself when your cell-level object lives somewhere else, for example the multi-omics union `preprocess/adata_sample.h5ad` (the multi-omics wrapper does exactly that) or your own file.

**Import:** `from sampledisco.sample_embedding.blocks import save_embedding_to_h5ad`

**Source:** `sample_embedding/blocks.py:283`

## Signature

```python
def save_embedding_to_h5ad(path: str, adata) -> str
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `path` | str | — | The `.h5ad` file. A missing file is written in full (returns `"rewrite"`). |
| `adata` | AnnData | — | The object holding both `uns['X_DR_sample']` and `uns['sample_embedding_params']`. |

## Returns

`"in-place"` when only the two `uns` entries were replaced, `"rewrite"` when the whole file was written.

!!! note "What the in-place check does not see"
    For `uns` entries other than the two embedding keys, only the key names are compared. If you changed the value of another `uns` entry in memory, the in-place path does not save that change.

## Output files

The file at `path`.

## Usage

Multi-omics chain, after [`cell_types_multiomics`](../multiomics/cell_types_multiomics_linux.md) wrote the union file. `compute_sample_embedding` cannot update this file on its own because it is not `preprocess/adata_preprocessed.h5ad`:

```python
import anndata as ad
from sampledisco.sample_embedding import compute_sample_embedding
from sampledisco.sample_embedding.blocks import save_embedding_to_h5ad

union_path = "api_demo/multiomics/preprocess/adata_sample.h5ad"
mo = ad.read_h5ad(union_path)
mo = compute_sample_embedding(mo, "api_demo/multiomics", modality_col="modality")
print(save_embedding_to_h5ad(union_path, mo))
```

It prints `in-place`: the file still matched the object, so only the two `uns` entries were replaced. The `compute_sample_embedding` call also rewrites `api_demo/multiomics/sample_embedding/`. About 7 s and 1.8 GB. On a path that does not exist it writes the whole object there (checked separately).
