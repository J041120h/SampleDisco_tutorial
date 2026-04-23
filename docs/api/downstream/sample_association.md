# `run_dimension_association_analysis`

Per-PC variance-explained decomposition of the sample embeddings against every metadata variable. For each sample embedding present in the pseudobulk object (`X_DR_expression`, `X_DR_proportion`) and every sample-level metadata column, the function fits a linear model of the component on the variable and records how much variance it explains. Continuous variables use design `[1, x]` (R² = squared Pearson correlation); categorical variables use `[1, one-hot(levels, drop-first)]` (R² = one-way ANOVA η²). Significance is assessed by permuting the variable and refitting, producing a directly-comparable metric across variable types and components.

Use this when you want to know which covariates (batch, age, study, sex, phenotype, ...) drive each PC of the sample embedding — it's the standard confounder / leading-covariate check before downstream tests.

**Source:** `sample_association/association.py:503`

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
| `pseudo_adata` | AnnData | — | Sample-level pseudobulk with `X_DR_expression` / `X_DR_proportion` in `.uns` or `.obsm` and per-sample metadata in `.obs`. |
| `output_dir` | str | — | Root directory. Created if missing. |
| `continuous_cols` | list, optional | `None` | Override automatic classification. Leave as `None` to auto-detect numeric columns. |
| `categorical_cols` | list, optional | `None` | Override automatic classification. Leave as `None` to auto-detect categorical columns. |
| `n_permutations` | int | `999` | Permutation count for the R² null. Set to `0` to skip. |
| `sample_col` | str | `"sample"` | Sample-id column excluded from testing. |
| `random_state` | int | `42` | Seed for the permutation RNG. |
| `verbose` | bool | `True` | Print progress and variable classification. |

## Returns

`dict` with keys:

| Key | Value |
| --- | --- |
| `"results"` | `{embedding_key: DataFrame}` of per-variable, per-component association tables. |
| `"continuous_cols"` | Final list of continuous variables tested. |
| `"categorical_cols"` | Final list of categorical variables tested. |

Each table has columns: `variable`, `component`, `kind`, `n_levels`, `r2`, `perm_p`, `n`, `pearson_r`, `spearman_r`, `fdr`.

## Output files

Under `{output_dir}/`:

- `variance_explained_expression.csv`, `variance_explained_proportion.csv` — the full R² tables.
- `figures/expression_variance_heatmap.pdf`, `figures/proportion_variance_heatmap.pdf` — component × variable R² heatmaps.
- `figures/expression_top_associations.pdf`, `figures/proportion_top_associations.pdf` — curated panels of the strongest associations.

## Usage

```python
from genodistance.sample_association import run_dimension_association_analysis

assoc = run_dimension_association_analysis(
    pseudo_adata=pseudo_adata,
    output_dir="/results/rna/sample_association",
    n_permutations=999,
    sample_col="sample",
    verbose=True,
)

expr_df = assoc["results"]["X_DR_expression"]
expr_df.query("perm_p < 0.05").sort_values("r2", ascending=False).head(10)
```
