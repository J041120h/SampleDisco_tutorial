# `transfer_labels`

k-nearest-neighbor label transfer in the sample embedding `X_DR_sample`. **References** are the samples with a label that are not in `query`; every other sample is a **query** (by default, the unlabelled samples). Each query takes the majority label of its `n_neighbors` nearest references (Euclidean distance on the first `n_pcs` PCs), optionally weighting votes by 1/distance; ties go to the label of the nearest tied neighbor. To estimate reliability, every reference is also labelled from the other references (leave-one-out; with `group_col`, references of the same group are left out too).

**Import:** `from sampledisco import transfer_labels` (or `from sampledisco.sample_prediction import transfer_labels`)

**Source:** `sample_prediction/prediction.py:137`

## Signature

```python
def transfer_labels(
    pseudo_adata: AnnData,
    label_col: str,
    *,
    query: Optional[Sequence[bool]] = None,
    n_neighbors: int = 5,
    weighted: bool = False,
    n_pcs: Optional[int] = None,
    group_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, Dict]
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `pseudo_adata` | AnnData | — | Sample-level AnnData ([`build_sample_adata`](../shared/build_sample_adata.md)) with `X_DR_sample` in `.uns` (DataFrame) or `.obsm`. |
| `label_col` | str | — | Label column in `.obs`. |
| `query` | bool array or None | `None` | Samples to label. `None`: the samples with a missing label. A mask can also hold out labelled samples. |
| `n_neighbors` | int | `5` | Number of nearest references that vote. |
| `weighted` | bool | `False` | Weight each vote by 1/distance. |
| `n_pcs` | int or None | `None` | Use the first `n_pcs` PCs; `None` = all. |
| `group_col` | str or None | `None` | Leave-one-out also drops references of the same group (a donor; `"sample"` for multi-omics units). |

## Returns

`(table, metrics)`:

- `table`: indexed by sample; `label` (input), `query`, `predicted` (the transferred label for queries, the leave-one-out label for references) and `confidence` (vote share of the assigned label).
- `metrics`: `loo_accuracy`, `loo_balanced_accuracy`, `n_reference`, `n_query`, `n_neighbors`, `n_pcs`.

## Usage

Continuing from the RNA chain; the first two samples are treated as unlabelled:

```python
import anndata as ad
import numpy as np
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco import transfer_labels

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")

query = np.zeros(pseudo.n_obs, dtype=bool)
query[:2] = True
table, metrics = transfer_labels(pseudo, "sev.level", query=query, n_neighbors=3)
print(table)
print(metrics)
```

Output (under a second):

```text
             label  query predicted  confidence
sample
CoV-3-Wilk       4   True         1    0.666667
CoV-34-Aruna     4   True         4    0.666667
CoV-36-Aruna     4  False         1    1.000000
CoV-4-Wilk       4  False         1    0.666667
HD-1-Wilk        1  False         1    0.666667
HD-20-Aruna      1  False         1    0.666667
HD-21-Aruna      1  False         1    0.666667
HD-4-Wilk        1  False         1    1.000000
{'loo_accuracy': 0.6666666666666666, 'loo_balanced_accuracy': 0.5, 'n_reference': 6, 'n_query': 2, 'n_neighbors': 3, 'n_pcs': 7}
```

`label` keeps the true value of the two held-out queries, so you can see that one of them (`CoV-3-Wilk`) is mislabelled. With only 2 severe references left, 3 neighbors mostly vote "1"; on 8 samples this is a demonstration of the call, not a test of the method.

See [Prediction and label transfer](../../tutorials/downstream/index.md#prediction) and its options there.
