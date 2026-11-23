# `run_autotune`

Chooses the **RMD weight α** of the sample embedding by optimizing a proxy score, instead of using the default `rmd_weight="equal"`. It builds the A1/A2/A3/RMD blocks of [`compute_sample_embedding`](calculate_sample_embedding.md) once, then for each candidate α re-weights the blocks, runs PCA and the sample-level batch correction, and scores the result. The best α is used to build the final embedding, which is written to `adata.uns['X_DR_sample']`. Only α is searched (`scope="alpha_only"`); the composition weights stay at their inverse-variance defaults.

The score uses two kinds of unit-level labels:

- **grouping**: `obs[grouping_col]` (first value per sample). Categorical → PC R² (one-way ANOVA on the PCs); numeric → CCA correlation.
- **batch**: the per-unit majority `batch_col` label (single-omics), the modality (multi-omics), or modality × batch (multi-omics with `batch_col`).

If there is no grouping column and only one batch label, nothing is searched: α falls back to the `"equal"` rule and the embedding is built once.

**Import:** `from sampledisco.parameter_selection.autotune import run_autotune`

**Source:** `parameter_selection/autotune.py:592`

## Signature

```python
def run_autotune(
    adata: AnnData,
    output_dir: str,
    *,
    sample_col: str = "sample",
    celltype_col: str = "cell_type",
    comp_emb_key: Optional[str] = None,
    rmd_emb_key: Optional[str] = None,
    modality_col: Optional[str] = None,
    batch_col: Optional[Union[str, List[str]]] = None,
    grouping_col: Optional[str] = None,
    medium_K: int = 120,
    fine_K: int = 300,
    rmd_dim: int = 8,
    pca_components: int = 10,
    batch_method: str = "harmony",
    scoring: str = "auto",
    search: str = "bayesian",
    scope: str = "alpha_only",
    alpha_bounds: Tuple[float, float] = (0.1, 100.0),
    seed: int = 42,
    save: bool = True,
    verbose: bool = True,
    tune_on_modality: Optional[str] = None,
    cluster_emb_key: Optional[str] = None,
    save_cell_adata: bool = True,
) -> Dict
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level object with cell types and both embeddings. Modified in place. |
| `output_dir` | str | — | Root output directory (see Output files). |
| `sample_col` / `celltype_col` | str | `"sample"` / `"cell_type"` | As in `compute_sample_embedding`. |
| `comp_emb_key` / `rmd_emb_key` | str or None | `None` | Embedding keys; `None` resolves them as in `compute_sample_embedding`. |
| `modality_col` | str or None | `None` | `"modality"` for multi-omics units. |
| `batch_col` | str, list or None | `None` | Sample-level batch column; only the first one is used (for the RMD group, the batch-mixing scores and the sample-level correction). |
| `grouping_col` | str or None | `None` | Phenotype the embedding should track (e.g. `"sev.level"`). Wrapper key: `*_autotune_grouping_col`. |
| `medium_K` / `fine_K` / `rmd_dim` / `pca_components` | int | `120` / `300` / `8` / `10` | Block and PCA settings, as in `compute_sample_embedding` (`rmd_dim` = `rmd_dim_per_cluster`). |
| `batch_method` | str | `"harmony"` | Sample-level correction used for every candidate and the final embedding. |
| `scoring` | str | `"auto"` | Proxy score, higher is better. See the table below. |
| `search` | str | `"bayesian"` | `"bayesian"` (Gaussian process with expected improvement: 5 log-spaced starts + 10 steps), `"golden_section"` (12 iterations on log α), or `"grid"` (15 log-spaced points). Any other value raises `ValueError`; there is no `"golden"`. |
| `scope` | str | `"alpha_only"` | The only supported value. |
| `alpha_bounds` | tuple | `(0.1, 100.0)` | Search range for α (searched on a log scale). |
| `seed` | int | `42` | Seed for blocks, PCA and Harmony. |
| `save` | bool | `True` | Write the CSV and the report. |
| `verbose` | bool | `True` | Print progress. |
| `tune_on_modality` | str or None | `None` | Multi-omics: score only the units of this modality (at least 5 needed) while the final embedding still covers all units. Output then goes to a separate folder (below). |
| `cluster_emb_key` | str or None | `None` | Deprecated spelling of `comp_emb_key`. |
| `save_cell_adata` | bool | `True` | Update `{output_dir}/preprocess/adata_preprocessed.h5ad` if it exists, as in `compute_sample_embedding`. |

**Scoring options.** Components that need a missing label (no grouping, or a single batch label) are left out of the ensembles; a single-metric scorer whose label is missing returns 0 for every α.

| `scoring` | Score |
| --- | --- |
| `"auto"`, `"multi_metric_proxy"` | Mean of min-max-scaled grouping tracking (PC R² or CCA), SPS (between/within-quartile distance ratio of the grouping), batch iLISI, and −ASW of batch. |
| `"ilisi_label"` | Mean of two terms: grouping tracking and batch iLISI. |
| `"cca"` | Grouping tracking only. |
| `"sps"`, `"cv_knn_severity"`, `"pseudotime_spearman"` | One grouping metric: SPS, negative 5-fold CV error of 3-NN regression on the grouping, \|Spearman(PC1, grouping)\|. |
| `"ilisi_batch"`, `"neg_asw_batch"` | One batch-mixing metric. |
| `"sev_minus_batch"` | Grouping tracking − 0.5 × batch ASW. |

## Returns

A `dict`:

| Key | Value |
| --- | --- |
| `best_params` | `{"rmd_weight": α}` |
| `best_score` | Score at α (`nan` when no search ran). |
| `trace` | List of `(α, score)` in evaluation order. |
| `block_weights` | `[w_A1, w_A2, w_A3, α]`. |
| `adata` | The input cell-level AnnData, with `uns['X_DR_sample']` and `uns['sample_embedding_params']` (here holding `best_params`, `best_score`, `search`, `scoring`, `K_c`, `K_med`, `K_fine`, `block_weights`, `n_evals`, `autotuned=True`, …). |
| `sample_embedding` | The final embedding as a `DataFrame` (units × PCs), the same as `uns['X_DR_sample']`. |

It does not return a sample-level AnnData; use [`build_sample_adata`](build_sample_adata.md).

## Output files

With `save=True`, under `{output_dir}/sample_embedding/` (or `{output_dir}/sample_embedding_tune-on-{modality}/` with `tune_on_modality`):

- `sample_embedding.csv`: the final embedding.
- `autotune_record.txt`: configuration, active proxies, best α and score, block weights, the α range within 1% of the best score, a note if α hit a bound, and the full search trace.

`{output_dir}/preprocess/adata_preprocessed.h5ad` is updated only if it exists and `save_cell_adata=True`. No `sample_embedding_blocks.npz` is written.

## Usage

Continuing from the RNA chain. A separate `output_dir` keeps the main run's files untouched; the in-memory `adata` gets the tuned embedding:

```python
import anndata as ad
from sampledisco.parameter_selection.autotune import run_autotune

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")
result = run_autotune(
    adata,
    "api_demo/rna_autotune",
    batch_col="batch",
    grouping_col="sev.level",
    search="grid",
)
print(result["best_params"], round(result["best_score"], 3))
print(result["sample_embedding"].shape)
```

**Writes** → `api_demo/rna_autotune/sample_embedding/sample_embedding.csv` and `autotune_record.txt`. Printed:

```text
{'rmd_weight': 100.0} 0.437
(8, 7)
```

All four `auto` proxies are active (grouping `sev.level`, batch `batch`). From `autotune_record.txt`:

```text
  best rmd_weight    : 100.000000
  best score        : 0.436552
  block weights     : 6.0828, 1.5706, 1.0000, 100.0000  (A1, A2, A3, RMD)
  total evaluations : 15
  α within 1% of best: [37.2759, 100]  (3/15 evaluations) — the objective is flat over this range
  NOTE: α equals a search bound; the optimum may lie outside [0.1, 100] — widen alpha_bounds and re-run.
```

The score goes from about 0.38 for α ≤ 1 to 0.437 at α = 100 and is flat above about 40. About 4 s and 2.0 GB.

!!! warning "Small cohorts"
    With 8 samples, every candidate embedding has at most 7 PCs and the proxies are flat or saturated, so the chosen α often lands on a bound of `alpha_bounds` (the report says so). Treat the result as a smoke test; see [Autotune options](../../tutorials/downstream/index.md#autotune).

!!! tip "Use the config-driven wrapper"
    Set `rna_autotune_enable` / `atac_autotune_enable` / `multiomics_autotune_enable: true` together with `*_derive_sample_embedding: true`; the wrapper then calls `run_autotune` **instead of** `compute_sample_embedding`. The keys `*_autotune_search`, `*_autotune_scoring`, `*_autotune_scope`, `*_autotune_alpha_bounds`, `*_autotune_grouping_col` and `multiomics_autotune_tune_on_modality` map onto the arguments above; `batch_col` comes from the sample-level batch key. See [Autotune](../../tutorials/downstream/index.md#autotune), including its options.
