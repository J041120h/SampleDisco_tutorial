# Prediction and label transfer

Two supervised uses of the sample embedding `uns['X_DR_sample']` (new in 0.5.0):

- `predict_phenotype` predicts a sample phenotype from the first *k* PCs with cross-validation, and reports how well it can be predicted.
- `transfer_labels` labels new or unlabelled samples from their nearest labelled neighbours in the embedding.

Both take the one-row-per-sample AnnData from `build_sample_adata` (see [Before you start](index.md#before-you-start)) and read labels from its `.obs`. They return tables and write no files. The config-driven wrapper does not run them; call them from Python.

## Phenotype prediction

```python
from sampledisco import predict_phenotype

pred, metrics = predict_phenotype(
    pseudo_adata,
    label_col="sev.level",   # binary or multiclass label in pseudo_adata.obs
    group_col=None,          # e.g. "donor" when a donor has several samples
    n_pcs=None,              # first k PCs; None = all
    task="classification",   # or "regression" for a continuous phenotype
    score_by="sample",       # or "group": score the per-group mean prediction
    n_splits=5,
    random_state=0,
)
print(metrics["auc"], metrics["balanced_accuracy"])
```

**Model.** The PCs are standardized. A classification label is fitted with an L2 logistic regression (binary or multiclass), a continuous one with ridge regression. The penalty is chosen by an inner 3-fold cross-validation inside each training split, scored by AUC (binary), balanced accuracy (multiclass) or R², from C in {10, 1, 0.1, 0.01} or alpha in {0.1, 1, 10, 100}; when several penalties score equally, the weakest is kept. The outer split is `n_splits`-fold, stratified by class and seeded by `random_state`. Each sample is predicted once, by the model trained without its fold. Samples with a missing label are skipped. This is the model used for the Sound Life trait prediction in the SampleDisco paper.

!!! warning "Use `group_col` when a donor has several samples"
    With repeated draws of a donor, or several units of one sample, pass the donor (or sample) column as `group_col`. Every sample of a group then stays in one fold, both in the outer and in the inner split. Without it, a model can recognise the donor rather than the phenotype, and the scores are inflated.

**Scoring level.** By default every sample is scored. With `score_by="group"`, the out-of-fold probabilities (or predictions) and the labels are first averaged within each `group_col` group (the majority label for classification), so each donor counts once. This is how the paper's donor-level AUCs were computed:

```python
pred, metrics = predict_phenotype(pseudo_adata, "cmv", group_col="donor", score_by="group")
print(metrics["auc"])   # one point per donor
```

**Returns** `(pred, metrics)`:

- `pred`: a DataFrame indexed by sample, with `label`, `predicted`, `fold` and, for classification, one `prob_<class>` column per class (out-of-fold). It is per sample for both `score_by` values.
- `metrics`: a dict with `score_by`, `n_samples`, `n_groups`, `n_pcs`, the penalty chosen in each fold (`C` or `alpha`), and `auc` (binary; the second class in sorted order is positive), `accuracy`, `balanced_accuracy`, or `r2`, computed at the `score_by` level.

## Label transfer

```python
from sampledisco import transfer_labels

table, metrics = transfer_labels(
    pseudo_adata,
    label_col="sev.level",   # samples with a missing value are the queries
    query=None,              # or a boolean mask: hold out labelled samples as queries
    n_neighbors=5,
    weighted=False,          # True: each neighbour votes with weight 1/distance
    n_pcs=None,
    group_col=None,          # e.g. "donor"; "sample" for multi-omics units
)
print(table[table["query"]])
print(metrics["loo_accuracy"])
```

Each query takes the majority label of its `n_neighbors` nearest labelled samples (Euclidean distance on the first `n_pcs` PCs). Ties go to the label of the nearest tied neighbour. The confidence is the share of the votes won by the assigned label (weighted by 1/distance with `weighted=True`).

To judge how far to trust the transfer, every labelled sample is also labelled from the other labelled samples (leave-one-out). With `group_col`, the other samples of the same group are left out as well. For multi-omics this matters: the RNA and ATAC units of one sample are usually each other's nearest neighbours, so pass `group_col="sample"`.

**Returns** `(table, metrics)`:

- `table`: a DataFrame indexed by sample, with `label` (the input), `query` (True for samples that were labelled by the transfer), `predicted` (the transferred label for queries, the leave-one-out label for labelled samples) and `confidence`.
- `metrics`: a dict with `loo_accuracy`, `loo_balanced_accuracy`, `n_reference`, `n_query`, `n_neighbors` and `n_pcs`.

!!! note "Small cohorts"
    Both scores are estimates from the labelled samples. With few samples per class (the 8-sample demo has 4 per class) they are very noisy; use them to check that the call works, not to judge the embedding.

See the API pages for [`predict_phenotype`](../../api/downstream/predict_phenotype.md) and [`transfer_labels`](../../api/downstream/transfer_labels.md).
