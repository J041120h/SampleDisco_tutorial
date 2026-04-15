# `CCA_Call`

Supervised pseudotime through Canonical Correlation Analysis. For each of the two sample embeddings (`X_DR_expression`, `X_DR_proportion`), CCA projects the embedding onto the one-dimensional axis most correlated with a phenotype column in `adata.obs` (e.g. severity, age, disease stage). The leading CCA axis is reported as a sample-level pseudotime and visualized as a 2D scatter using the PCs that contribute most to the axis. Useful when you already have a phenotype gradient and want to inspect how samples order along it.

**Source:** `sample_trajectory/CCA.py:388`

## Signature

```python
def CCA_Call(
    adata: AnnData,
    output_dir: str = None,
    trajectory_col: str = "sev.level",
    n_components: int = 2,
    auto_select_best_2pc: bool = True,
    verbose: bool = False,
    show_sample_labels: bool = False,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Sample-level AnnData with DR results in `.uns` (keyed by `X_DR_expression` / `X_DR_proportion`). |
| `output_dir` | str, optional | `None` | If provided, writes to `{output_dir}/CCA/`. |
| `trajectory_col` | str | `"sev.level"` | `.obs` column used as the CCA target. |
| `n_components` | int | `2` | PCA components drawn from DR before running CCA. Use `10` for RNA, typically `2` for ATAC. |
| `auto_select_best_2pc` | bool | `True` | Pick the best 2-PC combination for the 2D scatter. |
| `verbose` | bool | `False` | Print diagnostics. |
| `show_sample_labels` | bool | `False` | Overlay sample IDs on the plot (clutters when many samples). |

## Returns

`Tuple[float, float, Dict, Dict]` — `(proportion_score, expression_score, proportion_pseudotime, expression_pseudotime)`:

- `*_score` — CCA correlation coefficient along the leading axis.
- `*_pseudotime` — `{sample_id: float}` dict mapping each sample to its position along the CCA axis (used later by [`run_trajectory_gam_differential_gene_analysis`](trajectory_gam_dge.md)).

## Output files

Under `{output_dir}/CCA/`:

- `pca_{n_components}d_cca_expression.pdf`
- `pca_{n_components}d_cca_proportion.pdf`
- Per-DR pseudotime CSVs and PC-contribution diagnostics.

## Usage

```python
from genodistance.sample_trajectory import CCA_Call

prop_score, expr_score, prop_time, expr_time = CCA_Call(
    adata=pseudo_adata,
    output_dir="/results/rna",
    trajectory_col="sev.level",
    n_components=10,
    auto_select_best_2pc=True,
    verbose=True,
)
```
