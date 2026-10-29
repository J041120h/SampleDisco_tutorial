# `predict_phenotype`

Out-of-fold prediction of a sample phenotype from the first `n_pcs` PCs of the sample embedding `X_DR_sample` (0.5.0+). Standardized PCs feed an L2 logistic regression (binary or multiclass) or, for a continuous phenotype, a ridge regression. The penalty is chosen by an inner 3-fold grouped cross-validation inside each training split; the outer split is stratified, grouped and seeded. With `group_col`, every sample of a group (for example a donor) stays in one fold, so no donor is in training and test at once.

**Source:** `sampledisco/sample_prediction/prediction.py:31`

## Signature

```python
def predict_phenotype(
    pseudo_adata: AnnData,
    label_col: str,
    *,
    group_col: Optional[str] = None,
    n_pcs: Optional[int] = None,
    task: str = "classification",
    score_by: str = "sample",
    n_splits: int = 5,
    random_state: int = 0,
) -> Tuple[pd.DataFrame, Dict]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudo_adata` | AnnData | — | Sample-level AnnData (`build_sample_adata`) with `X_DR_sample` in `.uns` (DataFrame) or `.obsm`. |
| `label_col` | str | — | Phenotype column in `.obs`. Samples with a missing value are skipped. |
| `group_col` | str, optional | `None` | Column whose groups must never be split across folds (donor, or sample for multi-omics units). `None` = every sample is its own group. |
| `n_pcs` | int, optional | `None` | Use the first `n_pcs` PCs; `None` = all. |
| `task` | str | `"classification"` | `"classification"` (logistic regression, C in {10, 1, 0.1, 0.01}) or `"regression"` (ridge, alpha in {0.1, 1, 10, 100}). Ties go to the weakest penalty. |
| `score_by` | str | `"sample"` | `"group"` averages the out-of-fold probabilities (or predictions) and labels per `group_col` group before scoring (majority label for classification). |
| `n_splits` | int | `5` | Number of outer folds. |
| `random_state` | int | `0` | Seed of the outer split. |

## Returns

`Tuple[DataFrame, dict]` — `(pred, metrics)`.

- `pred`: indexed by sample; `label`, `predicted`, `fold`, and `prob_<class>` per class for classification (out-of-fold).
- `metrics`: `task`, `score_by`, `n_samples`, `n_groups`, `n_pcs`, the penalty chosen per fold (`C` or `alpha`), and `auc` (binary; second class in sorted order is positive), `accuracy`, `balanced_accuracy` (classification) or `r2` (regression).

## Usage

```python
from sampledisco import predict_phenotype

pred, metrics = predict_phenotype(pseudo_adata, "cmv", group_col="donor", score_by="group")
```

See the [tutorial](../../tutorials/downstream/prediction.md).
