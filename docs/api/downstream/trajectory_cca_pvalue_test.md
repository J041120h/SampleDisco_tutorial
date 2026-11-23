# `cca_pvalue_test`

Permutation test for a CCA correlation. It shuffles the phenotype `num_simulations` times, fits a one-component CCA between the **first `n_pcs` columns** of the embedding and the shuffled phenotype, and records `|corr|` each time. With k the number of permuted correlations at least as large as `input_correlation` and N = `num_simulations`, the p-value is (1 + k) / (N + 1): it is never 0, and its smallest possible value is 1 / (N + 1). (0.5.0 returned k / N, which could be 0.) The phenotype is converted to numbers and samples with a missing value are dropped, as in [`CCA_Call`](trajectory_cca_call.md). A histogram of the null and a text summary are written.

**Import:** `from sampledisco.sample_trajectory.CCA_test import cca_pvalue_test`

**Source:** `sample_trajectory/CCA_test.py:557`

## Signature

```python
def cca_pvalue_test(
    pseudo_adata: AnnData,
    column: str,
    input_correlation: float,
    output_directory: str,
    num_simulations: int = 1000,
    trajectory_col: str = "sev.level",
    n_pcs: int = 2,
    random_state: int = 0,
    verbose: bool = True,
) -> float
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudo_adata` | AnnData | — | Sample-level AnnData; `.uns[column]` holds the coordinates and `.obs` the phenotype. |
| `column` | str | — | `.uns` key of the embedding, normally `"X_DR_sample"`. |
| `input_correlation` | float | — | Observed correlation, e.g. the score from [`CCA_Call`](trajectory_cca_call.md). |
| `output_directory` | str | — | Output goes to `{output_directory}/CCA_test/`. |
| `num_simulations` | int | `1000` | Number of permutations. |
| `trajectory_col` | str | `"sev.level"` | Phenotype column. |
| `n_pcs` | int | `2` | Leading columns used for the null (capped at the number available). The wrapper passes its `*_n_cca_pcs`. |
| `random_state` | int | `0` | Seed of the permutations. |
| `verbose` | bool | `True` | Print warnings and the runtime. |

!!! tip "Pass the same number of PCs as `CCA_Call`"
    The null is a CCA on the first `n_pcs` columns, which is the same statistic `CCA_Call` reports with its default `auto_select_best_2pc=True`: a CCA on all its `n_components` PCs (the best PC pair is used only for its plot and pseudotime). Set `n_pcs` to the number of PCs `CCA_Call` actually used, i.e. `n_components` clamped to the number of embedding columns (at most units − 1). The wrapper does this (`n_pcs = *_n_cca_pcs`), and the example uses 2 on both sides. With `auto_select_best_2pc=False`, `CCA_Call` scores PC1 + PC2 only, so pass `n_pcs=2`.

## Returns

`float`: the permutation p-value (1 + k) / (N + 1). With N = 1000 the smallest possible value is 1/1001 ≈ 0.001.

!!! note "Degenerate null (0.6.0)"
    When the PCs used reach the number of samples − 1 (for example 7 PCs on 8 samples), CCA fits any labelling perfectly: the observed score is 1.0 by construction, every permuted correlation is 1.0 too, and the p-value is meaningless. 0.5.0 then crashed while drawing the histogram (`ValueError: Too many bins for data range`). 0.6.0 returns the p-value, warns (`RuntimeWarning: ... permutation null is degenerate ...; the p-value carries little information`) and draws a single-bin histogram. Keep `n_pcs` (and `CCA_Call`'s `n_components`) well below the number of samples.

## Output files

Under `{output_directory}/CCA_test/`:

- `cca_pvalue_distribution_{column}.png`: null histogram with the observed value.
- `cca_pvalue_result_{column}.txt`: observed correlation, p-value, the line `P-value definition: (1 + k) / (N + 1) with k = … of N = … permuted |r| >= observed`, and the number of samples used.

## Usage

Continuing from the RNA chain:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_trajectory.CCA import CCA_Call
from sampledisco.sample_trajectory.CCA_test import cca_pvalue_test

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

score, _, _, _ = CCA_Call(pseudo, trajectory_col="sev.level", n_components=2)
p = cca_pvalue_test(pseudo, "X_DR_sample", score, "api_demo/rna",
                    num_simulations=1000, trajectory_col="sev.level", n_pcs=2)
print(round(score, 3), p)
```

**Writes** → `api_demo/rna/CCA_test/cca_pvalue_distribution_X_DR_sample.png` and `cca_pvalue_result_X_DR_sample.txt`. Printed (about 2 s):

```text
[CCA] Flipped pseudotime so it increases with 'sev.level' (Spearman rho before flip: -0.873)
P-value for observed correlation 0.7738352008760679: 0.04195804195804196 [(1 + 41) / (1000 + 1)]
0.774 0.04195804195804196
```

`cca_pvalue_result_X_DR_sample.txt` holds the same numbers plus `P-value definition: (1 + k) / (N + 1) with k = 41 of N = 1000 permuted |r| >= observed` and `N samples used: 8`. 0.5.0 reported k / N = 0.041 for this run. With 8 samples and two severity levels, this p-value is fragile; see the [Sound Life walkthrough](../../tutorials/soundlife_rna.md) for a 32-sample test.
