# `CCA_Call`

Supervised pseudotime by canonical correlation analysis (CCA) between the sample embedding and a numeric phenotype in `.obs` (severity, age, stage, …):

1. Take the first `n_components` PCs of `X_DR_sample` (capped at the number available). `trajectory_col` is converted to numbers; samples whose value is not numeric (`NaN`) are dropped.
2. **Score.** With `auto_select_best_2pc=True` (the default), the score is the correlation of a one-component CCA on **all** `n_components` PCs. With `auto_select_best_2pc=False` it is the CCA correlation on PC1 + PC2 only.
3. **Plot and pseudotime.** A 2-PC pair is chosen for the scatter plot: with `auto_select_best_2pc=True` and more than 2 PCs, the pair with the highest CCA correlation among all C(n, 2) pairs; otherwise PC1 + PC2. The samples are projected on the CCA direction fitted in that pair, and the projection, min-max scaled to [0, 1], is the pseudotime. The sign of a CCA direction is arbitrary, so the pseudotime is flipped (t → 1 − t) when its Spearman correlation with `trajectory_col` is negative, and a line `[CCA] Flipped pseudotime so it increases with '<trajectory_col>' (Spearman rho before flip: …)` is printed. The pseudotime then increases with the phenotype. The score (an absolute correlation) is unchanged.

The best pair is therefore used only for the plot and the pseudotime, never for the score. Any error is caught: the function then returns `(nan, nan, {}, {})` (with `verbose=True` it prints the error).

**Import:** `from sampledisco.sample_trajectory.CCA import CCA_Call`

**Source:** `sample_trajectory/CCA.py:377`

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
| `adata` | AnnData | — | Sample-level AnnData with `X_DR_sample` in `.uns` (or `.obsm`) and `trajectory_col` in `.obs`. |
| `output_dir` | str or None | `None` | If set, plots and the pseudotime CSV go to `{output_dir}/CCA/`. |
| `trajectory_col` | str | `"sev.level"` | Phenotype column; converted with `pd.to_numeric(errors="coerce")`. |
| `n_components` | int | `2` | Number of leading PCs used, clamped to the number of embedding columns (at most units − 1). Wrapper key: `*_n_cca_pcs`. |
| `auto_select_best_2pc` | bool | `True` | `True`: score on all `n_components` PCs, and the best PC pair for the plot and pseudotime. `False`: PC1 + PC2 for both. Only matters when `n_components > 2`. |
| `verbose` | bool | `False` | Print progress and errors. |
| `show_sample_labels` | bool | `False` | Write sample ids on the scatter plot. |

## Returns

A 4-tuple `(score, score, pseudotime, pseudotime)`, kept for backward compatibility; both slots hold the same `X_DR_sample` result.

- `score` (float): the CCA correlation.
- `pseudotime` (dict): `{sample_id: value in [0, 1]}`, ready for [`run_trajectory_gam_differential_gene_analysis`](trajectory_gam_dge.md).

## Output files

Under `{output_dir}/CCA/`:

- `pca_{n_components}d_cca_sample.pdf`: the chosen PC pair with the CCA direction.
- `pca_{n_components}d_cca_sample_contributions.pdf`: PC contributions.
- `pseudotime_sample.csv`: columns `sample`, `pseudotime`.

## Usage

Continuing from the RNA chain. `sev.level` is 1 (healthy) or 4 (severe):

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_trajectory.CCA import CCA_Call

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

score, _, ptime, _ = CCA_Call(pseudo, output_dir="api_demo/rna", trajectory_col="sev.level", n_components=2)
print(round(score, 3))
print({s: round(float(t), 2) for s, t in ptime.items()})
```

**Writes** → `api_demo/rna/CCA/pca_2d_cca_sample.pdf`, `pca_2d_cca_sample_contributions.pdf` and `pseudotime_sample.csv`. Printed (about 1 s):

```text
[CCA] Flipped pseudotime so it increases with 'sev.level' (Spearman rho before flip: -0.873)
0.774
{'CoV-3-Wilk': 0.61, 'CoV-34-Aruna': 1.0, 'CoV-36-Aruna': 0.48, 'CoV-4-Wilk': 0.6, 'HD-1-Wilk': 0.29, 'HD-20-Aruna': 0.4, 'HD-21-Aruna': 0.31, 'HD-4-Wilk': 0.0}
```

The fitted CCA direction pointed from severe to healthy, so the pseudotime was flipped. All four severe samples now lie above all four healthy ones.

!!! warning "Small cohorts"
    The score is an in-sample CCA fit, so it grows with the number of PCs. Once the PCs used reach the number of samples − 1 (for example `n_components=10` on the 8-sample demo, clamped to its 7 PCs), CCA fits any labelling and the score is 1.0 by construction. Keep `n_components` well below the number of samples, and check the score with [`cca_pvalue_test`](trajectory_cca_pvalue_test.md), passing `n_pcs` equal to the number of PCs used here. See [Small cohorts](../../tutorials/downstream/index.md#small-cohorts).
