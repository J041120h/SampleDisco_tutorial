# `migrate_obsm_keys` and the key resolvers

SampleDisco stores two cell-level embeddings in `.obsm`: `Z_comp` (sample removed) and `Z_rmd` (sample preserved). Older releases used other names (`Z_clust`, `Z_cmd`, `X_pca_harmony`, …; full table on [Legacy embedding key names](../../legacy_keys.md)). The module `sampledisco.utils.embedding_keys` holds the alias tables and the functions every step uses to find the right key.

**Import:** `from sampledisco.utils.embedding_keys import migrate_obsm_keys, resolve_comp_key, resolve_rmd_key, resolve_embedding_keys`

**Source:** `utils/embedding_keys.py`: `resolve_comp_key` :96, `resolve_rmd_key` :116, `resolve_embedding_keys` :149, `migrate_obsm_keys` :175.

## `migrate_obsm_keys`

```python
def migrate_obsm_keys(adata, *, verbose: bool = True) -> List[Tuple[str, str]]
```

Copies each legacy embedding onto its canonical name, **in memory only**: `Z_clust`, `X_pca_harmony`, `X_lsi_harmony`, `X_glue_harmony` → `Z_comp`; `Z_cmd`, `X_pca_harmony_nosamp`, `X_lsi_harmony_nosamp` → `Z_rmd`. The first match wins for each target. Legacy keys are kept, an existing canonical key is never overwritten, and no file is written.

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Modified in place. |
| `verbose` | bool | `True` | Print one line per copied key. |

**Returns** the list of `(old, new)` pairs that were applied (empty if nothing changed).

## Resolvers

```python
def resolve_comp_key(adata, override=None, *, fallbacks=(), required=True, context="") -> Optional[str]
def resolve_rmd_key(adata, override=None, *, comp_key=None, fallbacks=(), required=True, context="") -> Optional[str]
def resolve_embedding_keys(adata, comp_override=None, rmd_override=None, *,
                           rmd_required=True, rmd_fallbacks=(), comp_fallbacks=(), context="") -> Tuple[str, Optional[str]]
```

| Argument | Meaning |
| --- | --- |
| `override` / `comp_override` / `rmd_override` | An explicit key. If it is in `.obsm` it is returned as is. If not: in 0.6.0, the canonical name `Z_comp` / `Z_rmd` falls back to the first legacy alias present, with a `FutureWarning`; any other missing key raises `KeyError`. |
| `fallbacks` / `comp_fallbacks` / `rmd_fallbacks` | Extra keys tried after the alias table (e.g. `("X_glue",)`). |
| `required` / `rmd_required` | `False` returns `None` instead of raising when nothing is found. |
| `comp_key` | RMD resolver only: when no sample-preserved key exists and `required=True`, return this key with a `RuntimeWarning` (the RMD block is then near-degenerate). |
| `context` | Tag used in warning messages. |

With no override, each resolver scans its alias table top to bottom (canonical name first; a legacy name gives a `FutureWarning`, once per key), then `fallbacks`. `resolve_embedding_keys` calls both and returns `(comp_key, rmd_key)`.

These are the rules behind every `None` default for `comp_emb_key`, `rmd_emb_key`, `cell_embedding_column`, `use_rep` and `embedding_key` in the API.

!!! warning "Explicit key names on old files"
    In 0.6.0, passing `comp_emb_key="Z_comp"` (or `use_rep="Z_comp"`, `rmd_emb_key="Z_rmd"`) to a file that only has the legacy names uses the legacy alias (e.g. `Z_clust`, `Z_cmd`) with a `FutureWarning`. On 0.5.0 the same call raises `KeyError`. Code meant for both versions should leave the keys at `None`, or call `migrate_obsm_keys(adata)` right after loading.

## Usage

The pre-integrated multi-omics demo file uses the pre-0.3 names:

```python
import anndata as ad
from sampledisco.utils.embedding_keys import migrate_obsm_keys, resolve_embedding_keys

mo = ad.read_h5ad("data/test_multiomics_integrated.h5ad")
print(sorted(mo.obsm.keys()))
print(resolve_embedding_keys(mo))      # legacy names, with a FutureWarning
print(migrate_obsm_keys(mo))
print(resolve_embedding_keys(mo))      # canonical names now
```

Output (after two `FutureWarning`s for `Z_clust` and `Z_cmd`):

```text
['X_glue', 'X_umap', 'Z_clust', 'Z_cmd']
('Z_clust', 'Z_cmd')
[sampledisco] migrated obsm['Z_clust'] → obsm['Z_comp']
[sampledisco] migrated obsm['Z_cmd'] → obsm['Z_rmd']
[('Z_clust', 'Z_comp'), ('Z_cmd', 'Z_rmd')]
('Z_comp', 'Z_rmd')
```

Nothing is written to disk; save the object yourself to keep the new names. Runs in about 2 s (0.7 GB).
