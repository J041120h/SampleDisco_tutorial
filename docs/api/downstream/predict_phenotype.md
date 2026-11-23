# `predict_phenotype`

Out-of-fold prediction of a sample phenotype from the first `n_pcs` PCs of `X_DR_sample` (0.5.0+). The standardized PCs feed an L2 logistic regression (binary or multiclass) or, with `task="regression"`, a ridge regression. In each outer training split, an inner 3-fold grouped cross-validation picks the penalty (scored by AUC for two classes, balanced accuracy for more, R² for regression; ties go to the weakest penalty). The outer split has `n_splits` folds and is seeded by `random_state`; for classification it is stratified by class. With `group_col` (e.g. a donor), all samples of a group stay in the same fold, outer and inner, so no donor is in training and test at once. Samples with a missing label are left out.

**Import:** `from sampledisco import predict_phenotype` (or `from sampledisco.sample_prediction import predict_phenotype`)

**Source:** `sample_prediction/prediction.py:31`

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
| `pseudo_adata` | AnnData | — | Sample-level AnnData ([`build_sample_adata`](../shared/build_sample_adata.md)) with `X_DR_sample` in `.uns` (DataFrame) or `.obsm`. |
| `label_col` | str | — | Phenotype column in `.obs`. |
| `group_col` | str or None | `None` | Column whose groups are never split across folds (a donor; `"sample"` for multi-omics units). `None`: each sample is its own group. |
| `n_pcs` | int or None | `None` | Use the first `n_pcs` PCs; `None` = all. |
| `task` | str | `"classification"` | `"classification"` (C in {10, 1, 0.1, 0.01}) or `"regression"` (alpha in {0.1, 1, 10, 100}). For regression the outer folds are grouped and shuffled but not stratified. |
| `score_by` | str | `"sample"` | `"group"`: average the out-of-fold probabilities (or predictions) and the labels over `group_col` before scoring (majority label for classification). |
| `n_splits` | int | `5` | Outer folds. For classification each class needs at least `n_splits` samples (groups), or scikit-learn raises `ValueError`. |
| `random_state` | int | `0` | Seed of the outer split. |

## Returns

`(pred, metrics)`:

- `pred`: indexed by sample; `label`, `predicted`, `fold`, and for classification one `prob_<class>` column per class (out-of-fold).
- `metrics`: `task`, `score_by`, `n_samples`, `n_groups`, `n_pcs`, the penalty chosen per fold (`C` or `alpha`), and `auc` (two classes; the second class in sorted order is positive), `accuracy`, `balanced_accuracy` (classification) or `r2` (regression).

## Usage

Continuing from the RNA chain. The demo has 4 samples per `sev.level` class, so `n_splits` is at most 4:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco import predict_phenotype

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

pred, metrics = predict_phenotype(pseudo, "sev.level", n_splits=4)
print(pred.round(2))
print({k: metrics[k] for k in ("auc", "accuracy", "balanced_accuracy", "C")})
```

Output (under a second):

```text
             label predicted  prob_1  prob_4  fold
sample
CoV-3-Wilk       4         1    0.54    0.46     0
CoV-34-Aruna     4         4    0.09    0.91     1
CoV-36-Aruna     4         4    0.05    0.95     2
CoV-4-Wilk       4         1    0.65    0.35     3
HD-1-Wilk        1         4    0.05    0.95     0
HD-20-Aruna      1         1    0.59    0.41     3
HD-21-Aruna      1         4    0.04    0.96     1
HD-4-Wilk        1         4    0.25    0.75     2
{'auc': 0.3125, 'accuracy': 0.375, 'balanced_accuracy': 0.375, 'C': [10, 10, 10, 10]}
```

With `sev.level` values 1 and 4, class `4` is the positive class for the AUC. Below chance here; an earlier run on another node gave AUC 0.5. Six training samples per fold are not enough to learn anything stable.

!!! warning "Eight samples"
    With 8 samples each test fold holds 2, so AUC and accuracy jump in steps and change with `random_state`. The [prediction tutorial](../../tutorials/downstream/prediction.md) and the [Sound Life walkthrough](../../tutorials/soundlife_rna.md) show it on 32 donors.
