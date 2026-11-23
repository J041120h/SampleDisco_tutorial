# Legacy embedding key names

SampleDisco writes two cell-level embeddings into `.obsm`, matching the two views in the paper's Methods:

| `.obsm` key | Paper symbol | Sample identity | Used for |
|---|---|---|---|
| `Z_comp` | $z_i^{\mathrm{comp}}$ | **removed** | cell typing, the A1/A2/A3 composition blocks, the latent anchors |
| `Z_rmd` | $z_i^{\mathrm{RMD}}$ | **preserved** | the RMD displacement block |

Objects written by older releases use other names. **You do not need to rewrite them**: every step that reads an embedding resolves the old names automatically.

## Alias table

Resolution goes top to bottom; the first key present in `.obsm` wins. A legacy name works with a `FutureWarning`. All legacy names will be removed in 1.0.

| View | Canonical key | Legacy names, in order |
|---|---|---|
| Composition (sample removed) | `Z_comp` | `Z_clust`, `X_pca_harmony` (RNA), `X_lsi_harmony` (ATAC), `X_glue_harmony` (multi-omics) |
| RMD (sample preserved) | `Z_rmd` | `Z_cmd`, `X_pca_harmony_nosamp` (RNA), `X_lsi_harmony_nosamp` (ATAC) |

Current releases also write a copy of `Z_comp` under `Z_clust`, so their files stay readable by older installs. The copy will be removed in 1.0.

## Explicit keys

`comp_emb_key` and `rmd_emb_key` default to `None`, which resolves through the table above; that is the recommended setting. With an explicit key:

- **The canonical name** (`"Z_comp"` or `"Z_rmd"`) falls back to its legacy alias, with a `FutureWarning`, when only the alias is present.
- **Any other key** that is not in `.obsm` raises `KeyError`.

Two keyword arguments were renamed; the old names still work with a `FutureWarning` until 1.0:

| Old | New |
|---|---|
| `cluster_emb_key=` | `comp_emb_key=` |
| `z_clust_key=` | `z_comp_key=` |

No config key changed: existing config files run unmodified.

## Migrating an object in memory

To normalize the keys on an object you have already loaded:

```python
from sampledisco.utils.embedding_keys import migrate_obsm_keys

applied = migrate_obsm_keys(adata)
print(applied)   # e.g. [('Z_clust', 'Z_comp'), ('Z_cmd', 'Z_rmd')]
```

This copies each legacy embedding onto its canonical name. The old keys stay, an existing canonical key is never overwritten, and no file is written: save the object yourself to keep the change. The demo file `test_multiomics_integrated.h5ad`, for example, has `Z_clust` and `Z_cmd`.

!!! warning "RMD blocks computed before 0.3.0 on `Z_cmd` objects"
    Through 0.2.0, an object with `Z_clust` + `Z_cmd` got its RMD block computed on the composition embedding, without a warning. Sample embeddings recomputed from such objects now differ; regenerate any metrics derived from them.
