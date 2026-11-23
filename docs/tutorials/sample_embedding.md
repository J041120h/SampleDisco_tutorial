# Sample embedding in depth

`compute_sample_embedding` turns a cell-level AnnData into one vector per sample. The RNA, ATAC, multi-omics and MaxFuse pipelines all end with it, and every downstream analysis reads its output. This page describes what it computes, what its options change, and runs it once on the COVID-19 demo.

```text
cell-level AnnData (obs: sample, cell_type, [batch], [modality])
 ├─ obsm['Z_comp']  sample removed ──┬─ A1  cell-type fractions             K_c columns
 │                                   ├─ A2  soft fractions, K_med anchors   K_med columns
 │                                   └─ A3  soft fractions, K_fine anchors  K_fine columns
 └─ obsm['Z_rmd']   sample preserved ── RMD per-cell-type displacement      ≤ K_c × rmd_dim_per_cluster columns
         │
         ├─ centre each block, scale it to weight w_A1, w_A2, w_A3, α
         ├─ concatenate → PCA (at most n_units − 1 PCs)
         └─ optional sample-level batch correction
         → uns['X_DR_sample']  (units × PCs)
```

In a config run, `*_derive_sample_embedding: true` calls it for you. Call it yourself, as at the end of this page, to try other options without redoing preprocessing: it takes seconds.

## Units

A unit is a row of the sample embedding.

- **Single-omics** (no `modality_col`): one unit per value of `sample_col`.
- **Multi-omics** (`modality_col="modality"`): one unit per (sample, modality) pair that has cells, named `<sample>_<modality>`, e.g. `HD-1-Wilk_RNA` and `HD-1-Wilk_ATAC`.

Each unit also belongs to a **group**, the pool of units the RMD block compares it with:

| Setting | Group of a unit |
|---|---|
| Multi-omics | its modality |
| Single-omics with `batch_col` | the batch most of the sample's cells carry (the first column if `batch_col` is a list) |
| Single-omics, `batch_col=None` (default) | one common group: every sample is compared with all the others |

The function needs at least 2 units and 2 cell types.

## Composition blocks A1, A2, A3

All three are computed on `Z_comp`, where Harmony has removed the sample, so a cell's position says *what state it is in*, not which sample it came from.

- **A1, cell-type composition:** the fraction of the unit's cells in each `celltype_col` label (`K_c` columns).
- **A2 and A3, anchor composition:** MiniBatchKMeans is fitted on the `Z_comp` of all cells pooled, with `K_med` anchors for A2 and `K_fine` for A3. Every cell is soft-assigned to every anchor with a Gaussian kernel (σ = the median cell-to-anchor distance), and a unit's row is the mean of its cells' assignments. These are compositions over data-driven states finer than the cell-type labels.
- **`use_clr=True`** replaces A1–A3 by their centred log-ratio. The default keeps plain fractions.

**Anchor caps.** The number of anchors is capped by the number of cells:

```text
K_med  = min(medium_K, max(2, N // 200))
K_fine = min(fine_K,   max(2, N // 100))
```

`N` is the **total number of cells passed in**, pooled over all samples and modalities, not the cells per sample. With the defaults `medium_K=120` and `fine_K=300`, the caps bind only below 24,000 cells (`K_med`) or 30,000 cells (`K_fine`). The COVID demo, with 29,988 cells after QC, gets `K_fine = 299`. The values used are stored in `uns['sample_embedding_params']`.

## The RMD block

RMD (reference-relative mean displacement) is computed on `Z_rmd`. This embedding keeps the differences between samples, so it shows *how* the cells of one type differ from sample to sample. For every unit `u` and cell type `k`:

- `d(u, k)` = the mean `Z_rmd` of u's type-k cells − the mean `Z_rmd` of all type-k cells in the **other units of u's group** (a leave-one-out reference);
- multiplied by `√n(u,k) / max_u √n(u,k)`, so displacements estimated from few cells are shrunk.

If no other unit of the group has type-k cells (for example, a sample alone in its batch), the reference is borrowed from the nearest other group. A unit without type-k cells gets `d = 0`. Each cell type's displacement matrix is then reduced by PCA to `min(rmd_dim_per_cluster, dim(Z_rmd), n_units − 1)` columns, and the pieces are concatenated. With 8 samples, each cell type contributes 7 columns.

## Block weights

Each block is centred and scaled to a Frobenius norm of `√n_units × w`, so its share of the total variance is proportional to `w²`. With `block_weights=None` (default):

| Block | Weight |
|---|---|
| A1 | `√(K_fine / K_c)` |
| A2 | `√(K_fine / K_med)` |
| A3 | `1` |
| RMD | `α`, set by `rmd_weight` |

- The coarse cell-type fractions weigh most: a composition block's variance is inversely proportional to its number of columns.
- **`rmd_weight="equal"`** (default) sets `α² = w_A1² + w_A2² + w_A3²`, so the RMD block carries as much variance as the three composition blocks together.
- **A number** fixes α (it must be positive). Lower values favour composition, higher values favour cell states. [Autotune](downstream/index.md#autotune) chooses α from the data.
- **`block_weights=[w_A1, w_A2, w_A3, w_RMD]`** overrides everything, including `rmd_weight`. It needs one weight per block (3 with `use_rmd=False`).

## PCA and sample-level batch correction

The stacked blocks are reduced by PCA to `min(pca_components, n_units − 1, number of features)` components, so the 8-sample demo gets 7 PCs, not 10.

Then, **only with at least 8 units and more than one batch label**, the PCs are batch-corrected. The labels are:

| `batch_col` | Single-omics | Multi-omics |
|---|---|---|
| `None` | none, so no correction | the modality |
| one column | the sample's majority batch | `<modality>__<batch>` |
| two or more columns | multi-covariate Harmony, one covariate per column | the same; list the modality if you want it |

`batch_method` must be one of:

- **`"harmony"`** (default): harmonypy, 30 iterations, seeded. If it fails, the function warns and falls back to `"linear"`.
- **`"linear"`**: regresses the one-hot batch labels out of every PC. It removes `n_batches − 1` directions, which can erase the biology when there are many batches.
- **`"none"`**: keeps the PCA as it is.

Any other value raises `ValueError`. Below about 20 units, Harmony has little or no effect: on the 8-sample demo, `"harmony"` and `"none"` give practically the same embedding.

## Which cell embeddings it reads

| Argument | `None` (default) resolves to, in order |
|---|---|
| `comp_emb_key` | `Z_comp`, then the legacy names `Z_clust`, `X_pca_harmony`, `X_lsi_harmony`, `X_glue_harmony` |
| `rmd_emb_key` | `Z_rmd`, then the legacy names `Z_cmd`, `X_pca_harmony_nosamp`, `X_lsi_harmony_nosamp` |

- The RNA and ATAC `preprocess` functions write both keys. In multi-omics, `harmonize_xglue` writes `Z_comp`, and `Z_rmd` is the scGLUE embedding `X_glue`.
- A legacy name works with a `FutureWarning` ([Legacy embedding key names](../api/utils/embedding_keys.md)).
- If no RMD embedding is found, the function falls back to the composition key with a `RuntimeWarning`. The RMD block is then near-degenerate, and the result is not the published method.

## Outputs

`compute_sample_embedding` modifies the cell-level AnnData in place and returns it, with two new entries:

- `uns['X_DR_sample']`: a DataFrame of units × PCs (`PC1`, `PC2`, …). Every downstream module reads this.
- `uns['sample_embedding_params']`: the values actually used: `K_c`, the capped `medium_K` and `fine_K`, `rmd_weight` (α) and `block_weights`, `pca_components` (the PCs used) and `pca_components_requested`, the embedding keys, `batch_col`, `batch_method`, `seed` and `backend` (`cpu` or `gpu`).

`build_sample_adata` turns these into a one-row-per-unit AnnData ([Outputs](outputs.md#load-the-results-back)). With `save=True` (default), the function also writes under `output_dir`:

| File | Content |
|---|---|
| `sample_embedding/sample_embedding.csv` | the same table as `uns['X_DR_sample']` |
| `sample_embedding/sample_embedding_blocks.npz` | the blocks before centring and weighting: `unit_ids`, `A1_cell_types`, `A1`, `A2`, `A3`, `RMD` |
| `preprocess/adata_preprocessed.h5ad` | updated with the two `uns` entries, **only if this file already exists** and `save_cell_adata=True` |

In a config run, the arguments come from these keys (RNA shown; ATAC and multi-omics are analogous):

| Argument | Config key |
|---|---|
| `medium_K`, `fine_K`, `rmd_dim_per_cluster` | `rna_sample_embedding_medium_K`, `rna_sample_embedding_fine_K`, `rna_sample_embedding_rmd_dim` |
| `use_clr`, `use_rmd`, `rmd_weight`, `block_weights` | `rna_sample_embedding_use_clr`, `…_use_rmd`, `…_rmd_weight`, `…_block_weights` |
| `pca_components`, `batch_method` | `rna_sample_embedding_pca_components`, `rna_sample_embedding_batch_method` |
| `batch_col` | `rna_sample_level_batch_col` (`atac_sample_level_batch_col`, `multiomics_batch_col`) |
| `save_cell_adata` | `rna_save_cell_adata_after_embedding` |

## Worked example: the COVID-19 demo

Start from the cell-level file of the [Quick start](quickstart.md) RNA run, which already has `Z_comp`, `Z_rmd` and 8 cell types, and recompute the embedding as the demo config does (`batch_col="batch"`, everything else default). A different `output_dir` leaves the Quick start files untouched:

```python
import anndata as ad
from sampledisco.sample_embedding import compute_sample_embedding

adata = ad.read_h5ad("sampledisco_demo_output/rna/preprocess/adata_preprocessed.h5ad")
adata = compute_sample_embedding(
    adata,
    output_dir="embedding_demo",
    sample_col="sample",
    celltype_col="cell_type",
    batch_col="batch",
)
```

```text
[sample_embedding] comp_emb=Z_comp, rmd_emb=Z_rmd
[sample_embedding] 8 units; 29988 cells; comp_emb dim=20
[A1] coarse cell-type composition: shape=(8, 8)
[A2/A3] MiniBatchKMeans K=120 and K=299...
[A2] shape=(8, 120)
[A3] shape=(8, 299)
[RMD] LOO displacement on rmd_emb...
  [RMD] shape=(8, 56)
[sample_embedding] weights=[6.114, 1.579, 1.0, 6.39] (n_blocks=4)
  [batch correction] 2 groups (group only)
[sample_embedding] wrote embedding_demo/sample_embedding/sample_embedding.csv
[sample_embedding] wrote embedding_demo/sample_embedding/sample_embedding_blocks.npz
[sample_embedding] done in 4.04s; shape=(8, 7)
```

harmonypy's own `INFO` lines are left out. What the log shows:

- **The anchor cap:** 29,988 cells give `K_fine = 299` instead of 300.
- **Small-sample limits:** each of the 8 cell types contributes 7 RMD columns (56 in all), and the embedding has 7 PCs.
- **The weights** are `[√(299/8), √(299/120), 1, α]`, with α = √(37.4 + 2.5 + 1) = 6.39 from `rmd_weight="equal"`.
- **Batch correction** used the two batches as labels; at 8 units it changes little.
- No h5ad was updated, because `embedding_demo/preprocess/adata_preprocessed.h5ad` does not exist.

## What the options do

| Option | Default | Effect | When to change it |
|---|---|---|---|
| `rmd_weight` | `"equal"` | Weight α of the RMD block. | Lower it when the phenotype is mainly a change in composition, raise it when it is a change in cell state. [Autotune](downstream/index.md#autotune) chooses it from the data. |
| `use_rmd` | `True` | `False` drops the RMD block: composition only. | To test whether a signal is pure composition. |
| `use_clr` | `False` | Centred log-ratio of A1–A3: changes in rare cell types weigh more. | To compare compositions on a log-ratio scale. |
| `medium_K`, `fine_K` | `120`, `300` | Number of anchors in A2 and A3, capped by the cell count. | Rarely; the A weights rescale automatically. |
| `batch_col` | `None` | Sample-level batch correction, **and** each sample's RMD reference becomes the other samples of its own batch. | A few large batches, each holding samples of every condition. |
| `batch_method` | `"harmony"` | `"linear"` regresses the batch out; `"none"` keeps the PCA. | Avoid `"linear"` with many batches. |
| `pca_components` | `10` | Number of PCs, capped at n_units − 1. | Only matters with more than 11 units. |

On the 8-sample demo the options make no measurable difference. Their effects show in larger cohorts, as below.

!!! note "A per-batch RMD reference can remove the donor signal"
    On the 32-donor Sound Life subset (19 processing batches, 11 of them with a single donor) we compared the options by how well the 10 PCs predict age group and CMV serostatus (cross-validated AUC). With the default `batch_col=None`, the AUC was 0.73 for age and 0.78 for CMV. With `batch_col="batch"`, the CMV AUC fell to 0.51, and it stayed at 0.52 with `batch_method="none"`. The loss came from the RMD reference, not from Harmony: each donor was compared only with its batch mates, and 11 donors had none. `batch_method="linear"` on 19 batches also erased age (AUC 0.29). Composition alone (`use_rmd=False`) predicted age best (0.88) and CMV not at all (0.48). When most batches hold one or a few samples, remove the batch at the cell level (`cell_level_batch_key`) and leave `batch_col=None`.
