# `summarize_h5ad`

Prints a structured overview of an `.h5ad` file: shape, layers, `obs`/`var` columns, `obsm`/`varm` keys, a check of `.X` (sparse or dense, dtype, density, integer-like or not, NaN/Inf, min/max), previews of `.X`, layers, `obs`, `var`, `obsm` and `varm`, and per-column summaries (dtype, missing values, unique count, top values or range). Use it before a run to check that `.X` holds raw counts and that the sample, batch and phenotype columns are where the config expects them (see [Using your own data](../../tutorials/own_data.md)).

The file is loaded fully into memory. Errors are caught and printed (`❌ Error reading …`) rather than raised.

**Import:** `from sampledisco.utils.inspector import summarize_h5ad`

**Source:** `utils/inspector.py:136`

## Signature

```python
def summarize_h5ad(
    h5ad_path: str,
    n_examples: int = 10,
    n_col_examples: int = 5,
    print_unique_column: str | None = None,
    unique_column_where: str = "auto",
    unique_max_values: int | None = 200,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `h5ad_path` | str | — | File to inspect. |
| `n_examples` | int | `10` | Rows shown in the `.X`, layer, `obs`, `var`, `obsm` and `varm` previews. |
| `n_col_examples` | int | `5` | Top values shown per `obs`/`var` column. |
| `print_unique_column` | str or None | `None` | Also list every unique value of this column. |
| `unique_column_where` | str | `"auto"` | Where to look for it: `"auto"` (obs, then var), `"obs"` or `"var"`. |
| `unique_max_values` | int or None | `200` | Refuse to list more unique values than this; `None` lists all. |

## Returns

`None`; everything is printed.

## Usage

```python
from sampledisco.utils.inspector import summarize_h5ad

summarize_h5ad("data/test_RNA.h5ad", n_examples=3, print_unique_column="sample")
```

The start of the output (about 2 s). The `.X` check confirms raw, integer-like counts:

```text
🔍 Loading AnnData from: data/test_RNA.h5ad

📦 Basic Info
  - Shape (cells × genes): 29989 × 14881
  - Layers: []
  - obs columns: ['celltype', 'sample', 'old', 'large', 'batch', 'sev.level']
  - var columns: []
  - obsm keys: []
  - varm keys: []

============================================================
🔎 Unique values for column 'sample' (found in adata.obs):
============================================================
  - n_unique: 8
     1. 'CoV-3-Wilk'
     ...
     8. 'HD-4-Wilk'

🔎 Inspecting .X matrix:
  - storage type: sparse (csc_matrix)
  - dtype: float64
  - .X nonzeros: 43474639 / 446266309 (9.7419% nonzero, 90.2581% zero)
  - value type: float — values appear integer-like (no fractional parts)
  - contains NaN: False
  - contains Inf: False
  - min value: 1
  - max value: 30539
```

It goes on with example values, cell names, and one block per `obs` column (dtype, missing values, unique count, top values).
