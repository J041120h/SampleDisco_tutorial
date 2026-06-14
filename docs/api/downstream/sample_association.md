# `run_dimension_association_analysis`

Per-PC variance-explained decomposition of the sample embedding against every metadata variable. For the single sample embedding stored in the pseudobulk object (`X_DR_sample`) and every sample-level metadata column, the function fits a linear model of the component on the variable and records how much variance it explains. Continuous variables use design `[1, x]` (R² = squared Pearson correlation); categorical variables use `[1, one-hot(levels, drop-first)]` (R² = one-way ANOVA η²). Significance is assessed by permuting the variable and refitting, producing a directly-comparable metric across variable types and components.

Use this when you want to know which covariates (batch, age, study, sex, phenotype, ...) drive each PC of the sample embedding — it's the standard confounder / leading-covariate check before downstream tests.

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
| `pseudo_adata` | AnnData | — | Sample-level pseudobulk with `X_DR_sample` in `.uns` (DataFrame) or `.obsm` and per-sample metadata in `.obs`. |
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
| `"results"` | `{embedding_key: DataFrame}` of per-variable, per-component association tables (single key `X_DR_sample`). |
| `"continuous_cols"` | Final list of continuous variables tested. |
| `"categorical_cols"` | Final list of categorical variables tested. |
| `"dropped_non_sample_level"` | Columns dropped because they vary within at least one sample (e.g. per-cell QC). |

Each table has columns: `variable`, `component`, `kind`, `n_levels`, `r2`, `perm_p`, `n`, `pearson_r`, `spearman_r`, `fdr`.

## Output files

Under `{output_dir}/`:

- `variance_explained_sample.csv` — the full R² table.
- `figures/sample_variance_heatmap.pdf` — component × variable R² heatmap.
- `figures/sample_top_associations.pdf` — curated panel of the strongest associations.

## Usage

```python
from sampledisco.sample_association.association import run_dimension_association_analysis

assoc = run_dimension_association_analysis(
    pseudo_adata=pseudo_adata,
    output_dir="/results/rna/sample_association",
    n_permutations=999,
    sample_col="sample",
    verbose=True,
)

sample_df = assoc["results"]["X_DR_sample"]
sample_df.query("perm_p < 0.05").sort_values("r2", ascending=False).head(10)
```
