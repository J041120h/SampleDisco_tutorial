# `run_dimension_association_analysis`

How much of each embedding PC is explained by each sample-level variable. For every (variable, PC) pair it fits a linear model of the PC on the variable and records R²: continuous variables use the design `[1, x]` (R² = squared Pearson correlation), categorical ones `[1, one-hot(levels, drop first)]` (R² = one-way ANOVA η²). Significance comes from permuting the variable and refitting, so continuous and categorical variables get directly comparable p-values; BH-FDR is applied across the table.

Use it to see which covariates (batch, age, sex, phenotype, …) drive which PCs before testing anything downstream.

**Import:** `from sampledisco.sample_association.association import run_dimension_association_analysis`

**Source:** `sample_association/association.py:541`

## Signature

```python
def run_dimension_association_analysis(
    pseudo_adata: AnnData,
    output_dir: str,
    continuous_cols: Optional[List[str]] = None,
    categorical_cols: Optional[List[str]] = None,
    n_permutations: int = 999,
    sample_col: str = "sample",
    random_state: int = 42,
    verbose: bool = True,
) -> dict
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudo_adata` | AnnData | — | Sample-level AnnData with `X_DR_sample` in `.uns` (DataFrame) or `.obsm`, and sample metadata in `.obs`. |
| `output_dir` | str | — | Output folder; created if missing. |
| `continuous_cols` | list or None | `None` | Continuous variables. `None`: detected automatically from `.obs`. |
| `categorical_cols` | list or None | `None` | Categorical variables. `None`: detected automatically. Each list is resolved on its own, so passing only one still auto-detects the other. |
| `n_permutations` | int | `999` | Permutations for the R² null; `0` skips them. |
| `sample_col` | str | `"sample"` | Sample column, excluded from testing. |
| `random_state` | int | `42` | Seed for the permutations. |
| `verbose` | bool | `True` | Print progress and the variable lists. |

Columns that vary within a sample are dropped. [`build_sample_adata`](../shared/build_sample_adata.md) already keeps only sample-level columns, so the automatic detection tests covariates such as `batch` and `sev.level`; pass explicit lists to choose which ones, and as continuous or categorical.

## Returns

`dict` with:

| Key | Value |
| --- | --- |
| `"results"` | `{"X_DR_sample": DataFrame}`, one row per (variable, PC) with columns `variable`, `component`, `kind`, `n_levels`, `r2`, `perm_p`, `n`, `pearson_r`, `spearman_r`, `fdr`. |
| `"continuous_cols"`, `"categorical_cols"` | The variables actually tested. |
| `"dropped_non_sample_level"` | Variables dropped because they vary within a sample. |

## Output files

Under `{output_dir}/`:

- `variance_explained_sample.csv`: the full table.
- `figures/sample_variance_heatmap.pdf`: PC × variable R² heatmap.
- `figures/sample_top_associations.pdf`: the strongest associations.

## Usage

Continuing from the RNA chain:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_association.association import run_dimension_association_analysis

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

assoc = run_dimension_association_analysis(
    pseudo,
    "api_demo/rna/sample_association",
    continuous_cols=["sev.level"],
    categorical_cols=["batch"],
    n_permutations=999,
)
df = assoc["results"]["X_DR_sample"]
print(df[["variable", "component", "r2", "perm_p", "fdr"]].sort_values("r2", ascending=False).head(6))
```

**Writes** → `api_demo/rna/sample_association/variance_explained_sample.csv` (14 rows: 2 variables × 7 PCs), `figures/sample_variance_heatmap.pdf` and `figures/sample_top_associations.pdf`. Printed (about 9 s, 2.1 GB):

```text
     variable component        r2  perm_p       fdr
13      batch       PC7  0.601255   0.022  0.196000
0   sev.level       PC1  0.577817   0.028  0.196000
4   sev.level       PC5  0.205550   0.383  0.945000
10      batch       PC4  0.197458   0.193  0.900667
3   sev.level       PC4  0.172404   0.299  0.945000
7       batch       PC1  0.108850   0.454  0.945000
```

With 8 samples nothing survives FDR; severity loads mostly on PC1, and the batch residual on the last PC.
