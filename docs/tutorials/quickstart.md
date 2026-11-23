# Quick start: the whole pipeline from one config

This page runs the complete SampleDisco pipeline on the COVID-19 demo data from a single YAML config. The data has 8 samples, 4 in each `sev.level` class. You run RNA alone first, from the command line and then from Python, and finish with RNA and ATAC together. Along the way you see what each stage does, where its results land, and how to load them back.

You need a working installation (see [Installation](../installation.md)). On 8 CPU cores the RNA-only run takes about 3 minutes and peaks at about 4 GB of RAM. The RNA + ATAC run takes 3 to 5.5 minutes, depending on the machine, and peaks at about 8.5 GB.

!!! note "Run everything from one folder"
    The demo config uses relative paths: it reads `./data/test_RNA.h5ad` and `./data/test_ATAC.h5ad` and writes to `./sampledisco_demo_output/`. Run every command and Python snippet on this page from the folder created in step 1.

## 1. Set up a run folder

```bash
mkdir -p sampledisco_quickstart/data
cd sampledisco_quickstart
curl -fL -o data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
sampledisco --init-config config.yaml
```

```text
Wrote the demo config to 'config.yaml'.
Edit the data paths / options, then run:
  sampledisco -m complex --config config.yaml
WARNING: 'config.yaml' is the config used for the Zenodo demo dataset — its QC thresholds and trajectory/grouping columns (e.g. 'sev.level', '*_supervised: true') are tuned for that demo data, not yours. Review and adjust every field before running on your own data.
```

- `test_RNA.h5ad` (94 MB) is all you need for steps 2–5. `test_ATAC.h5ad` (260 MB) is only used in step 6.
- Keep the `-f`: without it, curl saves an HTTP error page under the `.h5ad` name and still reports success. See [Demo data](demo_data.md) for other ways to get the files and for what they contain.
- `--init-config` writes the demo config that ships inside the package (the same settings as [`config_demo.yaml`](../assets/config_demo.yaml)). It refuses to overwrite an existing file. The warning is a reminder that the values are tuned for the demo data.

## 2. RNA-only run from the CLI

The template runs RNA and ATAC. For the first run, switch the ATAC pipeline off. That is the only edit.

```bash
sed -i.bak 's/^run_atac_pipeline: true/run_atac_pipeline: false/' config.yaml
grep -oE '^(output_dir|run_[a-z]+_pipeline|rna_cluster_dge|rna_cca_pvalue): [^ ]+' config.yaml
```

```text
output_dir: "./sampledisco_demo_output"
run_rna_pipeline: true
run_atac_pipeline: false
run_multiomics_pipeline: false
rna_cluster_dge: false
rna_cca_pvalue: false
```

`sed -i.bak` works with both GNU and macOS `sed` and keeps the original as `config.yaml.bak`. You can just as well change the line in a text editor. The `grep` also shows two stages the template leaves off, RAISIN (`rna_cluster_dge`) and the CCA p-value; both are explained below.

Now run the pipeline:

```bash
sampledisco -m complex --config config.yaml
```

- `-m complex` is the only mode and also the default, so `sampledisco --config config.yaml` does the same thing. The flag is kept for older scripts.
- Before running, the CLI checks the config against `wrapper()`: every key must be present and unknown keys are rejected. Only the three `*_save_cell_adata_after_embedding` keys may be missing. Change values; don't delete keys.
- The command exits with a non-zero status if a pipeline fails.
- On 8 CPU cores this run took 2 min 45 s and peaked at 4.0 GB of RAM.

A few lines of the log tell you most of what happened:

```text
Final shape: 29988 cells × 12301 genes
...
[cell_types] Found 20 clusters
  [cell_types] Aggregating 20 → 8 using dendrogram
...
[sample_embedding] 8 units; 29988 cells; comp_emb dim=20
...
  [batch correction] 2 groups (group only)
...
[CCA] Flipped pseudotime so it increases with 'sev.level' (Spearman rho before flip: -0.873)
...
CCA Analysis Summary:
  Sample: score=0.9883 (PC1 + PC3)
...
  → Successfully fitted 9991/15711 genes
...
  → Found 0 significant genes (FDR < 0.05)
...
RNA pipeline completed successfully!

Total execution time: 147.14 seconds (2.45 minutes)
```

Leiden clustering does not give the same number of clusters on every run (18–20 across our runs), so your counts may differ; they are always merged down to the 8 cell types the config asks for. The `2 groups` line means the sample-level Harmony step ran (see step 4). The `Flipped pseudotime` line means the CCA pseudotime was turned around so that it increases with `sev.level` (see step 3).

Writes → `sampledisco_demo_output/rna/`, one folder per stage (see step 4), and `sampledisco_demo_output/sys_log/main_process_status.json`.

!!! note "RAISIN is off because of memory"
    The template leaves RAISIN cluster DGE (`rna_cluster_dge`) off. On this demo, switching it on raised the RNA run's peak memory from 4 GB to about 15 GB: its ComBat step works on a dense genes × cells matrix. You can add it later without redoing the rest of the run; see [Add a stage later](outputs.md#add-a-stage-later-raisin).

!!! warning "Keep `rna_cca_pvalue: false` on the demo"
    The CCA score here is already close to 1 (0.988 in this run). With 8 samples, a CCA on 5 PCs lines up with a two-level label almost perfectly whatever the data, and with `rna_n_cca_pcs` at 7 or more (the number of samples minus one) the score is exactly 1. The permutation null is then degenerate and its p-value is noise; 0.5.0 could even crash while plotting it. From 0.6.0 the p-value is (1 + k) / (N + 1) for k of N permutations at least as extreme, so it is never reported as 0, but on this demo it still means nothing. Leave the p-value off on this demo. For a cohort large enough to test, see the 32-sample [Sound Life walkthrough](soundlife_rna.md), and see [CCA trajectory](downstream/trajectory_cca.md) for the settings.

## 3. The same run from Python

The CLI is a thin layer over one function, `wrapper()`, which takes the config keys as keyword arguments. The snippet below does the same run as step 2 and writes to the same `output_dir`, so you only need one of the two.

```python
import yaml
from sampledisco.wrapper.wrapper import wrapper

with open("config.yaml") as fh:
    cfg = yaml.safe_load(fh)
cfg["run_atac_pipeline"] = False    # the same edit as in step 2

res = wrapper(**cfg)
```

!!! note
    - Import the function from `sampledisco.wrapper.wrapper`; this works on 0.5.0 and every later release. The shorter `from sampledisco import wrapper` works from 0.6.0 on; on 0.5.0 it gives you a module that can't be called.
    - `wrapper()` does not check the key set the way the CLI does. A key you leave out silently takes the function's default, which is often not the demo value (for example `run_atac_pipeline=False`, `initialization=True`, `rna_min_cells=500`, `rna_proportion_test=False`). Start from the full template and change values in `cfg`.

The call took 2 min 50 s and peaked at 4.0 GB, about the same as the CLI run. `wrapper()` returns a dict. Look inside:

```python
print(sorted(res))
rna = res["rna_results"]
print(sorted(rna))
print("cell level:  ", rna["adata"].shape)
print("sample level:", rna["pseudo_adata"].shape)
print(rna["pseudo_adata"].obs[["sev.level", "batch", "pseudotime_sample"]])
print(res["status_flags"]["rna"])
```

```text
['rna_results', 'status_flags', 'system_info']
['adata', 'celltype_col', 'pseudo_adata', 'status_flags']
cell level:   (29988, 12301)
sample level: (8, 7)
              sev.level  batch  pseudotime_sample
sample                                           
CoV-3-Wilk            4   Wilk           0.619918
CoV-34-Aruna          4  Aruna           1.000000
CoV-36-Aruna          4  Aruna           0.791702
CoV-4-Wilk            4   Wilk           0.970249
HD-1-Wilk             1   Wilk           0.524531
HD-20-Aruna           1  Aruna           0.512458
HD-21-Aruna           1  Aruna           0.396701
HD-4-Wilk             1   Wilk           0.000000
{'preprocessing': True, 'cell_type_cluster': True, 'derive_sample_embedding': True, 'cca_based_cell_resolution_selection': False, 'sample_distance_calculation': True, 'trajectory_analysis': True, 'trajectory_dge': True, 'sample_cluster': True, 'proportion_test': True, 'cluster_dge': False, 'visualization': True, 'dimension_association_analysis': True}
```

| Key | What it holds |
| --- | --- |
| `status_flags` | For each pipeline (`rna`, `atac`, `multiomics`), which stages have finished, plus `system_info`. The same content is saved as `sys_log/main_process_status.json`. |
| `system_info` | Platform, Python version, `use_gpu`, and `gpu_available` (whether the GPU was actually used). |
| `rna_results` | `adata`: the cell-level AnnData after QC and cell typing, with the sample embedding in `uns['X_DR_sample']`. `pseudo_adata`: the sample-level AnnData, one row per sample, with the embedding as `X` and `obsm['X_DR_sample']`, and the per-sample metadata plus `pseudotime_sample` in `obs`. `celltype_col`: the cell-type column the embedding and the downstream stages used (`cell_type` once cell typing has run; new in 0.6.0). `status_flags`: the same dict as above. |
| `atac_results`, `multiomics_results` | The same three entries for the other pipelines, when they run. |
| `rna_error`, `atac_error`, `multiomics_error` | Present only with `continue_on_error: true`, holding the error message of a pipeline that failed. With the default `false`, the exception is raised instead. |

The embedding has 7 PCs, not the 10 that `rna_sample_embedding_pca_components` asks for: PCA on 8 samples gives at most 7. The pseudotime runs from 0 to 1 and increases with `rna_trajectory_col`, so the COVID samples (`sev.level` 4) sit at the high end. Since 0.6.0 the CCA pseudotime is flipped, with a `[CCA] Flipped pseudotime …` log line, whenever it would otherwise run against the trajectory variable; in 0.5.x its direction was arbitrary.

## 4. What each stage does

Each stage has an on/off flag in the config. The table lists them in the order they run, with the folder each writes under `sampledisco_demo_output/rna/`.

| Stage (flag) | What it does | Writes to |
| --- | --- | --- |
| `rna_preprocessing` | Filters genes and cells (`rna_min_cells`, `rna_min_genes`, `rna_pct_mito_cutoff`), normalizes, flags 2,000 HVGs and runs PCA. Harmony then runs twice: on `batch` + `sample`, giving `Z_comp` (sample differences removed), and on `batch` alone, giving `Z_rmd` (sample differences kept). | `preprocess/` |
| `rna_cell_type_cluster` | Leiden clustering on `Z_comp`. If it finds more clusters than `rna_n_target_cell_clusters` (8), they are merged along a dendrogram down to 8. Also computes a UMAP. | `preprocess/` |
| `rna_derive_sample_embedding` | The SampleDisco embedding: cell-type composition (A1), soft k-means composition at two resolutions (A2, A3), and the RMD block of within-cell-type shifts on `Z_rmd`. The blocks are weighted and stacked, then reduced by PCA. Harmony on `rna_sample_level_batch_col` (`batch`) then corrects the result at the sample level; it runs only with at least 8 samples in more than one batch. | `sample_embedding/` |
| `rna_sample_distance_calculation` | Cosine and correlation distances between samples, a heatmap, and a check of whether samples with the same `sev.level` sit closer together. | `Sample_distance/` |
| `rna_trajectory_analysis` | Supervised (`rna_trajectory_supervised: true`, as here): CCA between the first `rna_n_cca_pcs` (5) PCs and the numeric `rna_trajectory_col` (`sev.level`). The score uses all 5 PCs. The pseudotime is each sample's position, scaled to 0–1, along the CCA axis of the best-aligned pair of PCs, which is also the pair plotted. Unsupervised: TSCAN instead. | `CCA/` (`TSCAN/` when unsupervised) |
| `rna_trajectory_dge` | Runs only as part of the trajectory stage. Builds a per-sample pseudobulk and fits a spline GAM for each gene along the pseudotime. | `trajectoryDEG/` |
| `rna_sample_cluster` | K-means on the embedding, with `rna_cluster_number` (3) clusters. | `sample_cluster/` |
| `rna_proportion_test` | Tests whether cell-type proportions differ between groups (CLR transform + moderated t-test). The groups come from `rna_cluster_differential_gene_group_col`, which the demo sets to `sev.level`; only when that key is empty are the k-means clusters used. | `sample_cluster/proportion_test/` |
| `rna_cluster_dge` | RAISIN differential expression between the same groups, after ComBat on `batch`. Off in the template. | `raisin_results/` |
| `rna_visualize_data` | Cell-type dendrogram, and the sample embedding coloured by each of `rna_grouping_columns`. | `visualization/` |
| `rna_dimension_association_analysis` | For each embedding PC and each per-sample `obs` column: variance explained (R²), a permutation p-value (`rna_association_n_permutations`: 99) and FDR. | `sample_association/` |

The ATAC pipeline has the same stages under `atac_*` keys. Its preprocessing uses TF-IDF and LSI instead of normalization and PCA. [Understanding the outputs](outputs.md) lists every file each stage writes.

## 5. Where the results land

After step 2 or 3 the run folder looks like this:

```text
sampledisco_quickstart/
├── config.yaml
├── config.yaml.bak
├── data/
└── sampledisco_demo_output/
    ├── rna/
    │   ├── preprocess/           adata_preprocessed.h5ad, cell_type.csv, umap_cell_type.png
    │   ├── sample_embedding/     sample_embedding.csv, sample_embedding_blocks.npz
    │   ├── Sample_distance/      cosine/, correlation/
    │   ├── CCA/                  pseudotime_sample.csv, pca_5d_cca_sample*.pdf
    │   ├── trajectoryDEG/        gam_*_<timestamp>.tsv, visualizations/
    │   ├── sample_cluster/       kmeans_clusters_sample.csv, proportion_test/
    │   ├── sample_association/   variance_explained_sample.csv, figures/
    │   └── visualization/        cell_type_dendrogram.png, *_embedding_sev.level.png
    └── sys_log/
        └── main_process_status.json
```

The tables are plain CSV, so any tool can read them. For example, the first two PCs of the sample embedding next to the pseudotime and the k-means cluster of each sample:

```python
import pandas as pd

OUT = "sampledisco_demo_output/rna"
emb = pd.read_csv(f"{OUT}/sample_embedding/sample_embedding.csv", index_col=0)
ptime = pd.read_csv(f"{OUT}/CCA/pseudotime_sample.csv", index_col=0)
clusters = pd.read_csv(f"{OUT}/sample_cluster/kmeans_clusters_sample.csv", index_col=0)
print(emb[["PC1", "PC2"]].join(ptime).join(clusters).round(2))
```

```text
                PC1   PC2  pseudotime  cluster_sample_kmeans
sample                                                      
CoV-3-Wilk     4.44 -7.82        0.62                      0
CoV-34-Aruna -12.38 -5.74        1.00                      2
CoV-36-Aruna  -1.62  2.34        0.79                      0
CoV-4-Wilk    -5.90  5.39        0.97                      1
HD-1-Wilk      0.91  3.48        0.52                      0
HD-20-Aruna    1.57  1.65        0.51                      0
HD-21-Aruna    3.33  2.11        0.40                      0
HD-4-Wilk      9.66 -1.41        0.00                      0
```

![Sample embedding of the RNA run, coloured by sev.level](../resource/quickstart/rna_embedding_sev_level.png){ width="520" }

<div class="figure-caption"><code>rna/visualization/proportion_embedding_sev.level.png</code> from the run in step 3: the 8 samples on the first two embedding PCs, coloured by <code>sev.level</code> (red: COVID-19, 4; blue: healthy, 1).</div>

[Understanding the outputs](outputs.md) shows how to reload the cell-level AnnData and the sample-level AnnData, what is inside `sample_embedding_blocks.npz`, how to read the status file, and how to rerun only part of the pipeline.

## 6. RNA + ATAC in one run

`--init-config` won't overwrite `config.yaml`, so write the template again under a new name and run it unchanged. It runs RNA (again) and ATAC into the same `sampledisco_demo_output/`.

```bash
sampledisco --init-config config_rna_atac.yaml
sampledisco -m complex --config config_rna_atac.yaml
```

- On 8 CPU cores this took 5 min 21 s and peaked at 8.0 GB of RAM.
- The RNA results in `rna/` are overwritten by this run; ATAC results go to `atac/`, which has the same layout.
- The template turns ATAC trajectory DGE and RAISIN off (`atac_trajectory_dge`, `atac_cluster_dge`). Of the visualizations it draws only the cell-type dendrogram.
- ATAC sets no `atac_sample_level_batch_col`, so its embedding gets no sample-level Harmony; the log says `[batch correction] 1 groups`.
- The ATAC samples are named by SRA accession (`SRR14466464`, …), not by the RNA donor names, so the RNA and ATAC embeddings can't be matched sample by sample. To put both modalities in one embedding, use the [multi-omics pipeline](multiomics.md).

## 7. What 8 samples can and can't tell you

The demo is sized to run in minutes, not to support conclusions. Treat every statistic it produces as illustrative only:

- **Embedding:** at most 7 PCs.
- **Trajectory:** the supervised CCA fits whatever label you give it, and with 5 PCs for 8 samples the fit is close to perfect by construction. The pseudotime ordering is a picture, not a finding, and the permutation p-value is meaningless (see the warning in step 2).
- **Trajectory DGE:** the GAM could fit only about 10,000 of the ~15,700 genes (the rest are skipped for too few samples), and none is significant at FDR < 0.05. The package itself warns that a spline GAM is underpowered below about 12 samples.
- **Clusters:** k-means with 3 clusters on 8 samples leaves tiny clusters; in step 5, two of the three clusters hold a single sample.
- **Proportion test and RAISIN:** they compare 4 samples against 4.
- **Dimension association:** 99 permutations cannot give a p-value below 0.01, and across all the PC × column tests nothing survives FDR.
- **Phenotype prediction:** `predict_phenotype` is a separate Python call, not a pipeline stage. With 4 samples per class it needs `n_splits` ≤ 4, and scores from 8 samples are noisy (see [Prediction](downstream/prediction.md)).

For results that mean something, run the same config on a larger cohort: the [Sound Life walkthrough](soundlife_rna.md) uses 32 donors.

## 8. Next steps

- [Understanding the outputs](outputs.md): every output file, reloading results, and resuming a run.
- [Configuration guide](configuration.md): adapting the config to your own data.
- [The sample embedding in depth](sample_embedding.md): the blocks, their weights and the Harmony step, from Python.
- [Downstream analysis](downstream/index.md): each analysis as a standalone Python call.
- [Multi-omics pipeline](multiomics.md): RNA and ATAC in one joint embedding.
