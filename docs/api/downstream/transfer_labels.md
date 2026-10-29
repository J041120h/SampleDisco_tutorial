# `transfer_labels`

k-nearest-neighbour label transfer in the sample embedding `X_DR_sample` (0.5.0+). Query samples take the majority label of their `n_neighbors` nearest labelled samples (Euclidean distance on the first `n_pcs` PCs), optionally weighted by 1/distance; ties go to the label of the nearest tied neighbour. Every labelled sample is also labelled from the others (leave-one-out, or leave-group-out with `group_col`) to estimate how reliable the transfer is.

**Source:** `sampledisco/sample_prediction/prediction.py:137`

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
| `pseudo_adata` | AnnData | — | Sample-level AnnData (`build_sample_adata`) with `X_DR_sample` in `.uns` (DataFrame) or `.obsm`. |
| `label_col` | str | — | Label column in `.obs`. Labelled samples outside `query` are the references. |
| `query` | bool array, optional | `None` | Samples to label. `None` = the samples with a missing label; a mask can also hold out labelled samples. |
| `n_neighbors` | int | `5` | Number of nearest references that vote. |
| `weighted` | bool | `False` | Weight each vote by 1/distance. |
| `n_pcs` | int, optional | `None` | Use the first `n_pcs` PCs; `None` = all. |
| `group_col` | str, optional | `None` | In the leave-one-out estimate, also leave out references of the same group (donor; `"sample"` for multi-omics units). |

## Returns

`Tuple[DataFrame, dict]` — `(table, metrics)`.

- `table`: indexed by sample; `label` (input), `query`, `predicted` (transferred label for queries, leave-one-out label for references) and `confidence` (vote share of the assigned label).
- `metrics`: `loo_accuracy`, `loo_balanced_accuracy`, `n_reference`, `n_query`, `n_neighbors`, `n_pcs`.

## Usage

```python
from sampledisco import transfer_labels

table, metrics = transfer_labels(pseudo_adata, "sev.level", n_neighbors=5, group_col="donor")
```

See the [tutorial](../../tutorials/downstream/prediction.md).
