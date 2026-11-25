# `predict_phenotype`

Out-of-fold prediction of a sample phenotype from the first `n_pcs` PCs of `X_DR_sample`. Two models are available:

- **`model="mlp"`** (default): the standardized PCs feed a multilayer perceptron with one hidden layer of 100 ReLU units, trained with Adam (learning rate 0.001, batches of 200, 200 epochs) and an L2 weight penalty of 10. Nothing is tuned. With `task="regression"` the target is standardized too.
- **`model="knn"`**: each held-out sample takes the majority label (or, for regression, the mean) of its `n_neighbors` nearest training samples, by Euclidean distance on the raw PCs as in [`transfer_labels`](transfer_labels.md). The class probability is the neighbours' vote share. It fits no parameters, which suits small cohorts (tens of samples).

The outer split has `n_splits` folds and is seeded by `random_state`; for classification it is stratified by class. With `group_col` (e.g. a donor), all samples of a group stay in the same fold, so no donor is in training and test at once. Samples with a missing label are left out.

**Import:** `from sampledisco import predict_phenotype` (or `from sampledisco.sample_prediction import predict_phenotype`)

**Source:** `sample_prediction/prediction.py:53`

## Signature

```python
def predict_phenotype(
    pseudo_adata: AnnData,
    label_col: str,
    *,
    model: str = "mlp",
    group_col: Optional[str] = None,
    n_pcs: Optional[int] = None,
    task: str = "classification",
    n_neighbors: int = 5,
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
| `model` | str | `"mlp"` | `"mlp"` or `"knn"` (see above). |
| `group_col` | str or None | `None` | Column whose groups are never split across folds (a donor; `"sample"` for multi-omics units). `None`: each sample is its own group. |
| `n_pcs` | int or None | `None` | Use the first `n_pcs` PCs; `None` = all. |
| `task` | str | `"classification"` | `"classification"` (binary or multiclass) or `"regression"`. For regression the outer folds are grouped and shuffled but not stratified. |
| `n_neighbors` | int | `5` | Neighbours for `model="knn"`; capped at the size of each training split. |
| `score_by` | str | `"sample"` | `"group"`: average the out-of-fold probabilities (or predictions) and the labels over `group_col` before scoring (majority label for classification). |
| `n_splits` | int | `5` | Outer folds. For classification each class needs at least `n_splits` samples (groups), or scikit-learn raises `ValueError`. |
| `random_state` | int | `0` | Seed of the outer split and of the MLP. |

## Returns

`(pred, metrics)`:

- `pred`: indexed by sample; `label`, `predicted`, `fold`, and for classification one `prob_<class>` column per class (out-of-fold).
- `metrics`: `task`, `model`, `score_by`, `n_samples`, `n_groups`, `n_pcs`, `n_neighbors` (knn), and `auc` (two classes; the second class in sorted order is positive), `accuracy`, `balanced_accuracy` (classification) or `r2` (regression).

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
print({k: metrics[k] for k in ("model", "auc", "accuracy", "balanced_accuracy")})
```

Output (a few seconds):

```text
              label  predicted  prob_1  prob_4  fold
sample
CoV-3-Wilk        4          4    0.49    0.51     0
CoV-34-Aruna      4          4    0.48    0.52     1
CoV-36-Aruna      4          4    0.48    0.52     2
CoV-4-Wilk        4          4    0.49    0.51     3
HD-1-Wilk         1          4    0.49    0.51     0
HD-20-Aruna       1          4    0.49    0.51     3
HD-21-Aruna       1          4    0.45    0.55     1
HD-4-Wilk         1          4    0.49    0.51     2
{'model': 'mlp', 'auc': 0.625, 'accuracy': 0.5, 'balanced_accuracy': 0.5}
```

With `sev.level` values 1 and 4, class `4` is the positive class for the AUC. Trained on 6 samples per fold, the MLP's L2 penalty keeps every probability near 0.5, so all samples get the same label. For a cohort this small, the k-nearest-neighbour model is the alternative:

```python
pred, metrics = predict_phenotype(pseudo, "sev.level", model="knn", n_neighbors=3, n_splits=4)
print({k: metrics[k] for k in ("model", "n_neighbors", "auc", "accuracy", "balanced_accuracy")})
```

```text
{'model': 'knn', 'n_neighbors': 3, 'auc': 0.46875, 'accuracy': 0.625, 'balanced_accuracy': 0.625}
```

!!! warning "Eight samples"
    With 8 samples each test fold holds 2, so AUC and accuracy jump in steps and change with `random_state` and with the machine the embedding was built on. Neither score above means anything. See [Prediction and label transfer](../../tutorials/downstream/index.md#prediction); its options cover `group_col` and `score_by`.
