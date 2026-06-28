# `CCA_Call`

Supervised pseudotime through Canonical Correlation Analysis. CCA projects the single sample embedding (`adata.uns['X_DR_sample']`) onto the one-dimensional axis most correlated with a phenotype column in `adata.obs` (e.g. severity, age, disease stage). The leading CCA axis is reported as a sample-level pseudotime and visualized as a 2D scatter using the PCs that contribute most to the axis. Useful when you already have a phenotype gradient and want to inspect how samples order along it.

**Source:** `sample_trajectory/CCA.py:372`

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
| `adata` | AnnData | — | AnnData with the sample embedding in `.uns['X_DR_sample']` (falls back to `.obsm['X_DR_sample']`). |
| `output_dir` | str, optional | `None` | If provided, writes to `{output_dir}/CCA/`. |
| `trajectory_col` | str | `"sev.level"` | `.obs` column used as the CCA target. |
| `n_components` | int | `2` | PCA components drawn from DR before running CCA. Use `10` for RNA, typically `2` for ATAC. |
| `auto_select_best_2pc` | bool | `True` | Pick the best 2-PC combination for the 2D scatter. |
| `verbose` | bool | `False` | Print diagnostics. |
| `show_sample_labels` | bool | `False` | Overlay sample IDs on the plot (clutters when many samples). |

## Returns

`Tuple[float, float, Dict, Dict]` — a legacy 4-tuple `(score, score, pseudotime, pseudotime)`. With the single-key embedding, both score slots and both pseudotime slots hold the **same** `X_DR_sample` result; the 4-tuple shape is preserved for back-compat with the wrapper.

- `score` — CCA correlation coefficient along the leading axis.
- `pseudotime` — `{sample_id: float}` dict mapping each sample to its position along the CCA axis (used later by [`run_trajectory_gam_differential_gene_analysis`](trajectory_gam_dge.md)).

## Output files

Under `{output_dir}/CCA/`:

- `pca_{n_components}d_cca_sample.pdf` — 2D scatter with the CCA direction overlay.
- `pca_{n_components}d_cca_sample_contributions.pdf` — PC-contribution diagnostic plot.
- `pseudotime_sample.csv` — per-sample pseudotime (`sample`, `pseudotime` columns).

## Usage

```python
from sampledisco.sample_trajectory.CCA import CCA_Call

score, _, pseudotime, _ = CCA_Call(
    adata=adata,
    output_dir="sampledisco_demo_output/rna",
    trajectory_col="sev.level",
    n_components=10,
    auto_select_best_2pc=True,
    verbose=True,
)
```
