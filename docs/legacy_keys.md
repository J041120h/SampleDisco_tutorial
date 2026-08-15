# Legacy embedding key names

SampleDisco writes two cell-level embeddings into `.obsm`, matching the two views in the paper's Methods:

| `.obsm` key | Paper symbol | Sample identity | Used for |
|---|---|---|---|
| `Z_comp` | $z_i^{\mathrm{comp}}$ | **removed** | cell typing, the A1/A2/A3 composition blocks, the latent anchors |
| `Z_rmd` | $z_i^{\mathrm{RMD}}$ | **preserved** | the RMD displacement block |

Those are the canonical names as of **0.3.0**. Earlier releases used different names, and objects preprocessed by them are still on disk. **You do not need to rewrite them** — every read site in the package resolves the old names automatically.

## Alias table

Resolution goes top to bottom; the first key present in `.obsm` wins.

**Composition (sample-removed) view**

| Key | Written by | On read |
|---|---|---|
| `Z_comp` | 0.3.0+ | canonical, silent |
| `Z_clust` | 0.2.0 | works, `FutureWarning` |
| `X_pca_harmony` | pre-0.2 (RNA) | works, `FutureWarning` |
| `X_lsi_harmony` | pre-0.2 (ATAC) | works, `FutureWarning` |
| `X_glue_harmony` | pre-0.2 (multi-omics) | works, `FutureWarning` |

**RMD (sample-preserved) view**

| Key | Written by | On read |
|---|---|---|
| `Z_rmd` | 0.2.0+ | canonical, silent |
| `Z_cmd` | pre-0.2 | works, `FutureWarning` |
| `X_pca_harmony_nosamp` | pre-0.2 (RNA) | works, `FutureWarning` |
| `X_lsi_harmony_nosamp` | pre-0.2 (ATAC) | works, `FutureWarning` |

All legacy names are scheduled for removal in **1.0**.

!!! warning "0.2.0 silently mis-computed the RMD block on `Z_cmd` objects"
    Through 0.2.0 the RMD resolver checked only `Z_rmd`. On an object carrying `Z_clust` + `Z_cmd`, it fell through to the **composition** key, so the displacement block was computed on the sample-*removed* embedding — no error, no warning, near-degenerate block. 0.3.0 recognises `Z_cmd` and resolves it correctly, so **sample embeddings recomputed from those objects will differ from ones generated under 0.2.0.** If you have saved metrics derived from such an object, regenerate them.

## Forward compatibility

0.3.0 writes `Z_comp` **and** a duplicate under the old name `Z_clust`, so an h5ad produced by 0.3.0 is still readable by an already-installed 0.2.0. The duplicate is removed in 1.0.

## Migrating an object in memory

To normalize the keys on an object you have already loaded:

```python
from sampledisco.utils.embedding_keys import migrate_obsm_keys

applied = migrate_obsm_keys(adata)
print(applied)   # e.g. [('Z_clust', 'Z_comp'), ('Z_cmd', 'Z_rmd')]
```

This copies each legacy embedding onto its canonical name. It is **non-destructive** — the old keys stay in place, an existing canonical key is never overwritten, and **no file is written**. Save the object yourself if you want the change to persist.

## Renamed keyword arguments

| Old | New | Status |
|---|---|---|
| `cluster_emb_key=` | `comp_emb_key=` | accepted with a `FutureWarning`; removed in 1.0 |
| `z_clust_key=` | `z_comp_key=` | accepted with a `FutureWarning`; removed in 1.0 |

Both default to `None`, which means *auto-resolve through the table above* — the recommended setting. Passing an explicit key that is not in `.obsm` raises `KeyError` rather than silently falling back.

No YAML configuration key changed; existing config files run unmodified.
