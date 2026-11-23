# Prediction and label transfer

Two supervised uses of the sample embedding (new in 0.5.0):

- `predict_phenotype` predicts a sample phenotype from the first *k* PCs with cross-validation, and reports how well it can be predicted;
- `transfer_labels` labels new or unlabelled samples from their nearest labelled neighbours in the embedding.

Both take the one-row-per-sample AnnData from `build_sample_adata` and read the labels from its `.obs`. They return tables and write no files. The config-driven run does not call them; use them from Python.

The first part runs on the 8-sample COVID demo, to show the calls; run the [setup block](index.md#before-you-start) first. With 4 samples per class, the scores there are too noisy to mean anything. The second part runs on the 32-donor Sound Life subset, where they start to be informative.

## Phenotype prediction

```python
from sampledisco.sample_prediction import predict_phenotype, transfer_labels

pred, metrics = predict_phenotype(
    pseudo_adata,
    label_col="sev.level",    # binary or multiclass label in pseudo_adata.obs
    group_col=None,           # e.g. a donor column when a donor has several samples
    n_pcs=None,               # first k PCs; None = all
    task="classification",    # or "regression" for a continuous phenotype
    score_by="sample",        # or "group": score the per-group mean prediction
    n_splits=4,               # at most the size of the smallest class (4 here)
    random_state=0,
)
print(pred.round(2))
print({k: metrics[k] for k in ("n_samples", "n_pcs", "C", "auc", "balanced_accuracy")})
```

```text
              label  predicted  prob_1  prob_4  fold
sample
CoV-3-Wilk        4          4    0.09    0.91     0
CoV-34-Aruna      4          4    0.14    0.86     1
CoV-36-Aruna      4          4    0.03    0.97     2
CoV-4-Wilk        4          1    0.69    0.31     3
HD-1-Wilk         1          4    0.13    0.87     0
HD-20-Aruna       1          4    0.32    0.68     3
HD-21-Aruna       1          4    0.20    0.80     1
HD-4-Wilk         1          4    0.15    0.85     2
{'n_samples': 8, 'n_pcs': 7, 'C': [10, 10, 10, 10], 'auc': 0.6875, 'balanced_accuracy': 0.375}
```

Each model is trained on 6 samples and predicts 2. Here it calls almost every sample severe: the AUC of 0.69 says the severe samples tend to get higher probabilities, but the balanced accuracy of 0.38 is below chance. Nothing can be read from 8 samples.

`n_splits` must not exceed the number of samples in the smallest class; the default of 5 raises `ValueError: n_splits=5 cannot be greater than the number of members in each class` on the demo.

**Model.** The PCs are standardized. A classification label is fitted with an L2 logistic regression (binary or multiclass), a continuous one with ridge regression. The penalty is chosen by an inner 3-fold cross-validation inside each training split, scored by AUC (binary), balanced accuracy (multiclass) or R², from C in {10, 1, 0.1, 0.01} or alpha in {0.1, 1, 10, 100}; when several penalties score equally, the weakest is kept. The outer split is `n_splits`-fold, stratified by class for classification and seeded by `random_state`. Each sample is predicted once, by the model trained without its fold. Samples with a missing label are skipped. This is the model used for the Sound Life trait prediction in the SampleDisco paper.

!!! warning "Use `group_col` when a donor has several samples"
    With repeated draws of a donor, or several units of one sample, pass the donor (or sample) column as `group_col`. Every sample of a group then stays in one fold, both in the outer and in the inner split. Without it, a model can recognise the donor rather than the phenotype, and the scores are inflated.

**Scoring level.** By default every sample is scored. With `score_by="group"`, the out-of-fold probabilities (or predictions) and the labels are first averaged within each `group_col` group (the majority label for classification), so each donor counts once. This is how the paper's donor-level AUCs were computed.

**Returns** `(pred, metrics)`:

- `pred`: a DataFrame indexed by sample, with `label`, `predicted`, `fold` and, for classification, one `prob_<class>` column per class (out-of-fold). It is per sample for both `score_by` values.
- `metrics`: a dict with `task`, `score_by`, `n_samples`, `n_groups`, `n_pcs`, the penalty chosen in each fold (`C` or `alpha`), and `auc` (binary; the second class in sorted order is positive), `accuracy` and `balanced_accuracy`, or `r2`, computed at the `score_by` level.

## Label transfer

Here two labelled samples are hidden and treated as queries:

```python
query = pseudo_adata.obs_names.isin(["CoV-3-Wilk", "HD-4-Wilk"])   # pretend these are unlabelled
table, tmetrics = transfer_labels(
    pseudo_adata,
    label_col="sev.level",    # samples with a missing value are queries too
    query=query,
    n_neighbors=3,
    weighted=False,           # True: each neighbour votes with weight 1/distance
    n_pcs=None,
    group_col=None,           # e.g. a donor column; "sample" for multi-omics units
)
print(table)
print(tmetrics)
```

```text
              label  query  predicted  confidence
sample
CoV-3-Wilk        4   True          1    0.666667
CoV-34-Aruna      4  False          4    0.666667
CoV-36-Aruna      4  False          1    1.000000
CoV-4-Wilk        4  False          1    0.666667
HD-1-Wilk         1  False          1    0.666667
HD-20-Aruna       1  False          1    0.666667
HD-21-Aruna       1  False          1    0.666667
HD-4-Wilk         1   True          1    1.000000
{'loo_accuracy': 0.6666666666666666, 'loo_balanced_accuracy': 0.6666666666666666, 'n_reference': 6, 'n_query': 2, 'n_neighbors': 3, 'n_pcs': 7}
```

The two hidden samples are labelled from the 6 references: `HD-4-Wilk` correctly as healthy, `CoV-3-Wilk` wrongly as healthy (2 of its 3 nearest references are healthy). For the references, `predicted` is the leave-one-out label, correct for 4 of 6.

Each query takes the majority label of its `n_neighbors` nearest labelled samples (Euclidean distance on the first `n_pcs` PCs). Ties go to the label of the nearest tied neighbour. The confidence is the share of the votes won by the assigned label (weighted by 1/distance with `weighted=True`).

To judge how far to trust the transfer, every labelled sample is also labelled from the other labelled samples (leave-one-out). With `group_col`, the other samples of the same group are left out as well. For multi-omics this matters: the RNA and ATAC units of one sample are usually each other's nearest neighbours, so pass `group_col="sample"`.

**Returns** `(table, metrics)`:

- `table`: a DataFrame indexed by sample, with `label` (the input), `query` (True for samples that were labelled by the transfer), `predicted` (the transferred label for queries, the leave-one-out label for labelled samples) and `confidence`.
- `metrics`: a dict with `loo_accuracy`, `loo_balanced_accuracy`, `n_reference`, `n_query`, `n_neighbors` and `n_pcs`.

## A larger example: Sound Life

The [Sound Life](../demo_data.md) subset has 32 donors (16 young, 16 older; half CMV-seropositive) with 1,000 RNA cells each. Its sample metadata are `subject` (the donor), `age_group`, `age_at_draw`, `sex`, `cmv` and `batch`.

This block reuses the embedding from the [Sound Life RNA tutorial](../soundlife_rna.md) if you ran it in the same folder, and otherwise builds one: preprocessing with Harmony on `batch`, the AIFI cell types already in the file, and the sample embedding with `batch_col="batch"`. Building took 1 to 2 minutes and 4.3 GB of RAM on 8 CPUs; preprocessing and cell embedding are not bit-reproducible from run to run, so your scores below will differ a little from ours.

```python
import os
import anndata as ad
from sampledisco.preparation.rna_preprocess_cpu import preprocess
from sampledisco.sample_embedding import compute_sample_embedding
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

SL = "sampledisco_soundlife_output/rna"
path = f"{SL}/preprocess/adata_preprocessed.h5ad"
if not os.path.exists(path):
    preprocess(
        h5ad_path="data/soundlife_RNA.h5ad", sample_meta_path=None, output_dir=SL,
        sample_column="sample", cell_level_batch_key=["batch"], sample_level_batch_key="batch",
        min_cells=100, min_genes=200, pct_mito_cutoff=20, num_cell_hvgs=2000,
        cell_embedding_num_PCs=20, num_harmony_iterations=20,
    )
sl = ad.read_h5ad(path)
if "X_DR_sample" not in sl.uns:
    sl = compute_sample_embedding(
        sl, SL, sample_col="sample", celltype_col="cell_type",   # AIFI L2 labels
        batch_col="batch", use_gpu=False, save=True,
    )
sl_pseudo = build_sample_adata(sl, sample_col="sample")
print(sl_pseudo.shape, sl_pseudo.obs["subject"].nunique(), "donors")
```

```text
(32, 10) 32 donors
```

Predict CMV serostatus and age group with donor-level grouping and scoring, and age as a number:

```python
for label in ["cmv", "age_group"]:
    _, m = predict_phenotype(sl_pseudo, label, group_col="subject", score_by="group",
                             n_splits=5, random_state=0)
    print(f"{label:10s} AUC {m['auc']:.2f}  balanced accuracy {m['balanced_accuracy']:.2f}  "
          f"(n = {m['n_groups']} donors)")
_, m = predict_phenotype(sl_pseudo, "age_at_draw", task="regression", group_col="subject",
                         n_splits=5, random_state=0)
print(f"age_at_draw  R2 {m['r2']:.2f}")
```

```text
cmv        AUC 0.59  balanced accuracy 0.59  (n = 32 donors)
age_group  AUC 0.80  balanced accuracy 0.72  (n = 32 donors)
age_at_draw  R2 0.16
```

Age group is predictable from this RNA-only embedding (AUC 0.80 on held-out donors); CMV serostatus is not (AUC 0.59, close to chance), and age as a number only weakly (R² 0.16). Over `random_state` 0 to 9 on the same embedding, the age-group AUC ranged from 0.73 to 0.87, the CMV AUC from 0.30 to 0.59, and the age R² from 0.16 to 0.42. A rebuilt embedding moves the numbers too: in an earlier run of this page, the same block gave 0.83 for age group and 0.55 for CMV.

In this subset every donor has exactly one sample, so `group_col="subject"` and `score_by="group"` give the same numbers as the defaults. They matter on the full Sound Life cohort, where most donors have several visits: there, leaving them out lets a model recognise donors across folds.

!!! note "How noisy are these scores?"
    Each is an estimate from 32 out-of-fold predictions. As the ranges above show, the fold split alone can move an AUC by 0.1 or more. Try several values of `random_state` before reading much into a difference between two embeddings.

See the API pages for [`predict_phenotype`](../../api/downstream/predict_phenotype.md) and [`transfer_labels`](../../api/downstream/transfer_labels.md).
