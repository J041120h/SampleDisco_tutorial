# Autotune (RMD weight α)

The sample embedding stacks four feature blocks before its PCA: three cell-type **composition** blocks built on `Z_comp` (cell types, then about 120 and about 300 k-means states) and one **RMD displacement** block built on `Z_rmd` (how each sample's cells sit inside their cell type). The three composition weights are fixed by the block sizes. The one free knob is **α**, the weight of the RMD block (`rmd_weight`). The default, `"equal"`, gives the RMD block the same energy as the three composition blocks together.

`run_autotune` searches α and nothing else: `scope="alpha_only"` is the only scope it accepts. For each candidate α it rebuilds the embedding from blocks computed once, scores it, and keeps the best. It is optional; most runs keep `"equal"`.

Run the [setup block](index.md#before-you-start) first.

## What it scores

The score needs at least one of two sample-level columns:

- **`grouping_col`** (biology to keep), for example `sev.level`. Its proxies: how well the embedding tracks the column (the canonical correlation for a numeric column, the mean R² of the top PCs on the one-hot labels for a categorical one), and SPS, a between/within-group distance ratio.
- **`batch_col`** (technical signal to remove). Its proxies: iLISI, the k-nearest-neighbour mixing of batches, and the negative silhouette width of the batch labels.

With neither column, there is nothing to score: autotune logs "no search performed", uses α = `"equal"` and still writes the outputs.

| `scoring` | What is maximized |
| --- | --- |
| `"auto"` (default) = `"multi_metric_proxy"` | mean of the available proxies above, each min-max scaled |
| `"ilisi_label"` | two terms only: grouping tracking and batch iLISI |
| `"cca"`, `"sps"`, `"cv_knn_severity"`, `"pseudotime_spearman"` | one grouping proxy (`cv_knn_severity`: cross-validated kNN error on the grouping; `pseudotime_spearman`: absolute Spearman correlation of PC1 with the grouping) |
| `"ilisi_batch"`, `"neg_asw_batch"` | one batch proxy |
| `"sev_minus_batch"` | grouping tracking − 0.5 × batch silhouette |

| `search` | Evaluations over `alpha_bounds` (searched on a log10 scale) |
| --- | --- |
| `"bayesian"` (default) | 5 log-spaced points including both bounds, then 10 picked by expected improvement |
| `"golden_section"` | 14 interior points, assuming one peak |
| `"grid"` | 15 log-spaced points including both bounds |

!!! warning "α is chosen against `grouping_col`"
    The search rewards separation of `grouping_col`. Judging the tuned embedding with the same column is circular; check it on a label that took no part in the search.

## Call

`run_autotune` replaces `adata.uns['X_DR_sample']`, writes `<output_dir>/sample_embedding/sample_embedding.csv`, and re-saves `<output_dir>/preprocess/adata_preprocessed.h5ad` when that file exists. With `output_dir=OUT`, every later page would use the tuned embedding. To keep the default embedding while you look, pass a copy and a separate folder:

```python
from sampledisco.parameter_selection.autotune import run_autotune

tuned = run_autotune(
    adata.copy(),                        # the copy keeps adata's embedding unchanged
    output_dir="sampledisco_demo_output/rna_autotune",
    sample_col="sample",
    celltype_col="cell_type",
    batch_col="batch",                   # batch proxies
    grouping_col="sev.level",            # grouping proxies
    scoring="auto",
    search="bayesian",
    scope="alpha_only",
    alpha_bounds=(0.1, 100.0),
    seed=42,
)
print(sorted(tuned))
print(tuned["best_params"], round(tuned["best_score"], 3))
print("block weights (A1, A2, A3, RMD):", [round(w, 3) for w in tuned["block_weights"]])
for alpha, score in sorted(tuned["trace"]):
    print(f"  alpha = {alpha:8.3f}   score = {score:.3f}")
```

```text
[autotune.build_blocks] K-means K=120...
[autotune.build_blocks] K-means K=299...
[autotune] search=bayesian  scoring=auto  scope=alpha_only  bounds=(0.1, 100.0)
[autotune] best rmd_weight=100.0000  score=0.4437  (15 evals)
  [batch correction] 2 groups (group only)
[autotune] wrote sampledisco_demo_output/rna_autotune/sample_embedding/sample_embedding.csv
[autotune] wrote sampledisco_demo_output/rna_autotune/sample_embedding/autotune_record.txt
[autotune] done in 5.00s
['adata', 'best_params', 'best_score', 'block_weights', 'sample_embedding', 'trace']
{'rmd_weight': 100.0} 0.444
block weights (A1, A2, A3, RMD): [6.114, 1.579, 1.0, 100.0]
  alpha =    0.100   score = 0.383
  alpha =    0.192   score = 0.383
  alpha =    0.223   score = 0.383
  alpha =    0.562   score = 0.383
  alpha =    2.073   score = 0.390
  alpha =    3.162   score = 0.396
  alpha =    7.723   score = 0.418
  alpha =   12.362   score = 0.428
  alpha =   17.783   score = 0.433
  alpha =   20.983   score = 0.434
  alpha =   37.652   score = 0.437
  alpha =   50.049   score = 0.437
  alpha =   76.873   score = 0.438
  alpha =   84.502   score = 0.444
  alpha =  100.000   score = 0.444
```

The score rises from 0.38, where the composition blocks dominate (α below 1), to about 0.44 once the RMD block dominates (α above about 40), and is nearly flat from there. For scale, the default `"equal"` corresponds to α = 6.39 on this data (the square root of the sum of the squared composition weights). The search ended on its upper bound, α = 100. The last steps are tiny and not stable: an identical call in another job returned α = 84.5 instead, scoring α = 100 at 0.438. On 8 samples, read the trace as "the score prefers an RMD-dominated embedding, by a small margin", not as an exact α.

`run_autotune` returns a dict:

| Key | Content |
| --- | --- |
| `best_params` | `{"rmd_weight": α}` |
| `best_score` | the score at that α (`nan` when no search ran) |
| `trace` | every evaluated `(α, score)` pair, in order |
| `block_weights` | the four block weights at the chosen α (A1, A2, A3, RMD) |
| `adata` | the AnnData you passed, now with the tuned `uns['X_DR_sample']` and `uns['sample_embedding_params']` |
| `sample_embedding` | the tuned embedding as a DataFrame (samples × PCs) |

To use the tuned embedding downstream, build the sample-level object from it:

```python
pseudo_tuned = build_sample_adata(tuned["adata"], sample_col="sample")
print(pseudo_tuned.uns["sample_embedding_params"]["best_params"])
```

```text
{'rmd_weight': 100.0}
```

### When α lands on a bound

The `"bayesian"` and `"grid"` searches evaluate both ends of `alpha_bounds`. When the reported α equals one of them, the best value may lie outside the range: widen `alpha_bounds` and search again. `autotune_record.txt` then adds a line starting with `NOTE`, and when several α score within 1% of the best it adds their range (`α within 1% of best`). That is what happened above, so we searched one decade higher:

```python
with open("sampledisco_demo_output/rna_autotune/sample_embedding/autotune_record.txt") as fh:
    print("".join(line for line in fh if line.lstrip().startswith(("best", "NOTE", "α within"))))

wide = run_autotune(
    adata.copy(),
    output_dir="sampledisco_demo_output/rna_autotune_wide",
    sample_col="sample",
    celltype_col="cell_type",
    batch_col="batch",
    grouping_col="sev.level",
    alpha_bounds=(1.0, 1000.0),          # moved up by a factor of 10
)
print(wide["best_params"], round(wide["best_score"], 3))
print([(round(a), round(s, 3)) for a, s in sorted(wide["trace"])])
```

```text
  best rmd_weight    : 100.000000
  best score        : 0.443699
  α within 1% of best: [84.502, 100]  (2/15 evaluations) — the objective is flat over this range
  NOTE: α equals a search bound; the optimum may lie outside [0.1, 100] — widen alpha_bounds and re-run.

[autotune.build_blocks] K-means K=120...
[autotune.build_blocks] K-means K=299...
[autotune] search=bayesian  scoring=auto  scope=alpha_only  bounds=(1.0, 1000.0)
[autotune] best rmd_weight=1000.0000  score=0.4382  (15 evals)
  [batch correction] 2 groups (group only)
[autotune] wrote sampledisco_demo_output/rna_autotune_wide/sample_embedding/sample_embedding.csv
[autotune] wrote sampledisco_demo_output/rna_autotune_wide/sample_embedding/autotune_record.txt
[autotune] done in 4.74s
{'rmd_weight': 1000.0} 0.438
[(1, 0.384), (2, 0.391), (6, 0.41), (15, 0.431), (21, 0.434), (32, 0.436), (71, 0.438), (111, 0.438), (178, 0.438), (210, 0.438), (288, 0.438), (377, 0.438), (461, 0.438), (727, 0.438), (1000, 0.438)]
```

From α ≈ 70 to 1000 the score does not move: the RMD block already dominates the embedding, and more weight changes nothing. The wider search lands on its new upper bound only because all those points tie to three decimals. When widening moves the answer to the new bound **and** the score stays flat, there is no optimum to find: choose any α on the plateau, or keep the default. When the score keeps rising, widen again.

**Writes** → `sampledisco_demo_output/rna_autotune/sample_embedding/`:

- `sample_embedding.csv`: the tuned embedding;
- `autotune_record.txt`: settings, active proxies, the result and the full search trace.

With `tune_on_modality` (multi-omics only) the folder is `sample_embedding_tune-on-<modality>/` instead. Nothing else is written: no pseudobulk or block files.

## From the config

Autotune runs inside the sample-embedding stage, **instead of** the plain `compute_sample_embedding` call. It needs three things in the config:

- `rna_derive_sample_embedding: true`: with the stage off, autotune never runs;
- `rna_autotune_enable: true`;
- a grouping column or a sample-level batch column to score: `rna_autotune_grouping_col`, or `rna_sample_level_batch_col`. With neither, there is no search. With only the batch column, which is the demo config's default state, only the batch proxies are scored: α is then chosen for batch mixing alone.

Starting from the demo config (with `run_atac_pipeline: false`, as in the [Quick start](../quickstart.md)), set:

```yaml
rna_derive_sample_embedding: true
rna_autotune_enable: true
rna_autotune_grouping_col: "sev.level"
rna_sample_level_batch_col: "batch"
rna_autotune_search: "bayesian"
rna_autotune_scoring: "auto"
rna_autotune_scope: "alpha_only"
rna_autotune_alpha_bounds: [0.1, 100.0]
```

and run the pipeline:

```bash
sampledisco --config config.yaml
```

The tuned embedding replaces the default one in `rna/sample_embedding/sample_embedding.csv` and in `rna/preprocess/adata_preprocessed.h5ad`, `autotune_record.txt` lands next to it, and every downstream stage uses it. In our config run all four proxies were active (`has batch column: True`, `has grouping col: True`). The search returned α = 100, the upper bound (score 0.429), and its record showed both extra lines:

```text
  α within 1% of best: [17.7828, 100]  (7/15 evaluations) — the objective is flat over this range
  NOTE: α equals a search bound; the optimum may lie outside [0.1, 100] — widen alpha_bounds and re-run.
```

It is the same plateau as above, on an embedding from a fresh preprocessing run; an earlier config run returned α = 65.1, and another one put the 1% range at 37.7 to 100.

The ATAC and multi-omics branches have the same keys with `atac_` and `multiomics_`; multi-omics adds `multiomics_autotune_tune_on_modality`, which is passed as `tune_on_modality`. With autotune on, `*_sample_embedding_rmd_weight`, `*_sample_embedding_block_weights`, `*_sample_embedding_use_clr` and `*_sample_embedding_use_rmd` are ignored, and the embedding is built on the CPU even with `use_gpu: true`.

Each Python search (15 evaluations) took 3 to 10 s on 8 CPUs, and the session peaked at 2.8 GB with two copies of `adata` in memory. The config run took 4 to 6.5 minutes and 3.0 GB in total, for all RNA stages.

See the [`run_autotune` API page](../../api/shared/run_autotune.md) for every parameter.
