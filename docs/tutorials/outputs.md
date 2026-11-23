# Understanding the outputs

This page picks up where the [Quick start](quickstart.md) ends. It maps every stage of a config run to the files it writes, shows how to load the results back into Python, and explains how to rerun only part of the pipeline. Run the snippets from the `sampledisco_quickstart` folder, after the quick start's RNA run (steps 2 or 3) and its RNA + ATAC run (step 6).

Everything lands under `output_dir` (`./sampledisco_demo_output` in the demo config). There is one folder per pipeline: `rna/`, `atac/` and `multiomics/`, which you can move elsewhere with `rna_output_dir`, `atac_output_dir` and `multiomics_output_dir`. There is also a `sys_log/` folder.

## The output folder at a glance

This is the tree after following the quick start and then every recipe on this page, trimmed to one entry where a folder holds several files of the same kind. Sizes are from that run.

```text
sampledisco_demo_output/
├── sys_log/
│   └── main_process_status.json                  1.7 KB
├── rna/
│   ├── preprocess/
│   │   ├── adata_preprocessed.h5ad               685 MB
│   │   ├── cell_type.csv                         731 KB
│   │   └── umap_cell_type.png                    1.9 MB
│   ├── sample_embedding/
│   │   ├── sample_embedding.csv
│   │   └── sample_embedding_blocks.npz
│   ├── Sample_distance/
│   │   ├── cosine/
│   │   │   ├── distance_matrix_sample_DR.csv
│   │   │   ├── distance_check_results_sample_DR_cosine.txt
│   │   │   ├── distance_statistics_summary_cosine.csv
│   │   │   ├── sample_DR_coordinates.csv
│   │   │   └── sample_distance_sample_DR_heatmap.pdf
│   │   └── correlation/                          the same five files
│   ├── CCA/
│   │   ├── pseudotime_sample.csv
│   │   ├── pca_5d_cca_sample.pdf
│   │   └── pca_5d_cca_sample_contributions.pdf
│   ├── trajectoryDEG/
│   │   ├── gam_all_genes_<timestamp>.tsv         800 KB; one set of four per run (6 sets here)
│   │   ├── gam_significant_<timestamp>.tsv
│   │   ├── gam_pseudoDEGs_<timestamp>.tsv
│   │   ├── gam_summary_<timestamp>.txt
│   │   ├── differential_gene_result.txt
│   │   └── visualizations/                       01_results_summary.png, 04_volcano_plot.png,
│   │                                             05_gene_curves.png, 06_sample_density.png,
│   │                                             08_sample_level_curves.png
│   ├── sample_cluster/
│   │   ├── kmeans_clusters_sample.csv
│   │   ├── kmeans_sample_embedding.png
│   │   └── proportion_test/
│   │       ├── proportion_test_1_vs_4.csv
│   │       ├── proportion_test_significant_summary.txt
│   │       └── proportion_*.png                  4 heatmaps
│   ├── raisin_results/                           only after "Add a stage later: RAISIN"
│   │   ├── 1_vs_4/                               raisin_results.csv (1.0 MB), volcano.png, volcano_labeled.png
│   │   └── summary_plots/                        all_results_combined.csv, raisin_summary.csv/.txt, 3 PNGs
│   ├── sample_association/
│   │   ├── variance_explained_sample.csv
│   │   └── figures/                              sample_variance_heatmap.pdf, sample_top_associations.pdf
│   └── visualization/
│       ├── cell_type_dendrogram.png
│       ├── proportion_embedding_sev.level.png
│       └── expression_embedding_sev.level.png
└── atac/                                         same layout as rna/, but no trajectoryDEG/ or raisin_results/
    ├── preprocess/adata_preprocessed.h5ad        1.9 GB
    └── visualization/cell_type_dendrogram.png    the only ATAC plot the demo config asks for
```

The two `adata_preprocessed.h5ad` files take up almost all of the 2.5 GB. Everything else is small tables and figures.

## Stage by stage

Each stage writes into its own folder under `rna/` (or `atac/`). A stage that is switched off writes nothing, and its folder from an earlier run stays as it was.

| Stage (flag) | Folder | Files |
| --- | --- | --- |
| Preprocessing (`*_preprocessing`) and cell typing (`*_cell_type_cluster`) | `preprocess/` | `adata_preprocessed.h5ad`: the cell-level AnnData (below), rewritten by each of the two stages and updated in place by the embedding stage. `cell_type.csv`: `cell_id`, `cell_type`. `umap_cell_type.png`. |
| Sample embedding (`*_derive_sample_embedding`) | `sample_embedding/` | `sample_embedding.csv`: samples × PCs. `sample_embedding_blocks.npz`: the four blocks before weighting. Also stores `uns['X_DR_sample']` and `uns['sample_embedding_params']` in `adata_preprocessed.h5ad`, unless `*_save_cell_adata_after_embedding: false`. |
| Sample distance (`*_sample_distance_calculation`) | `Sample_distance/<method>/` | `distance_matrix_sample_DR.csv`, `sample_distance_sample_DR_heatmap.pdf`, `sample_DR_coordinates.csv`, `distance_statistics_summary_<method>.csv`, and `distance_check_results_sample_DR_<method>.txt` (mean distance within vs. between `sev.level` groups). |
| Trajectory (`*_trajectory_analysis`) | `CCA/` | `pseudotime_sample.csv` (`sample`, `pseudotime`), `pca_5d_cca_sample.pdf` and `pca_5d_cca_sample_contributions.pdf`, where 5 is `rna_n_cca_pcs`. An unsupervised run writes TSCAN results to `TSCAN/` instead. |
| Trajectory DGE (`*_trajectory_dge`) | `trajectoryDEG/` | `gam_all_genes_<timestamp>.tsv`, `gam_significant_<timestamp>.tsv`, `gam_pseudoDEGs_<timestamp>.tsv`, `gam_summary_<timestamp>.txt`, `differential_gene_result.txt`, and `visualizations/0*_*.png`. Each run adds a new set of timestamped files. |
| Sample clustering (`*_sample_cluster`) | `sample_cluster/` | `kmeans_clusters_sample.csv` (`sample`, `cluster_sample_kmeans`), `kmeans_sample_embedding.png`. |
| Proportion test (`*_proportion_test`) | `sample_cluster/proportion_test/` | `proportion_test_<a>_vs_<b>.csv` for each pair of groups (here `1_vs_4`, the two `sev.level` values), `proportion_test_significant_summary.txt` and four heatmaps. |
| RAISIN (`*_cluster_dge`) | `raisin_results/` | `<a>_vs_<b>/raisin_results.csv` plus volcano plots, and `summary_plots/` with the combined results, a summary and heatmaps. |
| Visualization (`*_visualize_data`) | `visualization/` | `cell_type_dendrogram.png` (Ward linkage of the per-cell-type mean `Z_comp`). `proportion_embedding_<col>.png` and `expression_embedding_<col>.png` for each grouping column: both are the same scatter of the sample embedding's first two PCs (the names are historical). |
| Dimension association (`*_dimension_association_analysis`) | `sample_association/` | `variance_explained_sample.csv`, one row per PC × per-sample `obs` column; `figures/sample_variance_heatmap.pdf`, `figures/sample_top_associations.pdf`. |
| Every run | `sys_log/` | `main_process_status.json` (see [The status file](#the-status-file)). |

The distance stage also writes a working file, `summary_sample.csv`, into the pipeline folder and deletes it at the end of the run. If `rna_summary_sample_csv_path` names a file that already exists, the run updates that file and, from 0.6.0 on, keeps it (0.5.0 deleted it). A file the run creates at that path is still deleted.

## Reload the cell-level AnnData

`preprocess/adata_preprocessed.h5ad` holds the cell-level data after QC, with the cell embeddings, the cell types and the sample embedding. Start Python in the run folder:

```python
import anndata as ad

OUT = "sampledisco_demo_output/rna"
adata = ad.read_h5ad(f"{OUT}/preprocess/adata_preprocessed.h5ad")
print(adata)
```

```text
AnnData object with n_obs × n_vars = 29988 × 12301
    obs: 'celltype', 'sample', 'old', 'large', 'batch', 'sev.level', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'cell_type_original'
    var: 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'highly_variable_clust', 'highly_variable_rmd'
    uns: 'X_DR_sample', 'X_name', 'cell_type', 'cell_type_colors', 'cluster_mapping', 'hvg', 'log1p', 'neighbors', 'pca', 'sample_embedding_params', 'umap'
    obsm: 'X_pca', 'X_pca_rmd', 'X_umap', 'Z_clust', 'Z_comp', 'Z_rmd'
    varm: 'PCs'
    layers: 'counts'
    obsp: 'connectivities', 'distances'
```

| Field | Content |
| --- | --- |
| `X` | log-normalized expression of all genes that passed QC (12,301 here) |
| `layers['counts']` | the raw counts |
| `var['highly_variable']` | the 2,000 HVGs behind `X_pca` and `Z_comp`; `highly_variable_rmd` flags the set behind `X_pca_rmd` and `Z_rmd` |
| `obs['cell_type']` | the 8 cell types from the cell-typing stage. `celltype` is the label that came with the demo file, and `cell_type_original` holds the Leiden clusters before they were merged |
| `obsm['Z_comp']` | cell embedding with sample differences removed (Harmony on `batch` + `sample`); feeds cell typing and the composition blocks. `Z_clust` is a copy under the 0.2.0 name, kept so older versions can still read the file ([legacy keys](../legacy_keys.md)) |
| `obsm['Z_rmd']` | cell embedding with sample differences kept (Harmony on `batch` only); feeds the RMD block |
| `uns['X_DR_sample']` | the sample embedding, a DataFrame of samples × PCs |
| `uns['sample_embedding_params']` | the settings that produced it (see below) |

The file is about 700 MB on the demo, because it keeps every gene twice (`X` and `counts`).

## Rebuild the sample-level AnnData

The downstream analyses work on a small AnnData with one row per sample. The pipeline builds it in memory from the cell-level file; `build_sample_adata` does the same for you:

```python
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

pseudo_adata = build_sample_adata(adata, sample_col="sample")
print(pseudo_adata)
print(pseudo_adata.obs[["sev.level", "batch"]])
```

```text
AnnData object with n_obs × n_vars = 8 × 7
    obs: 'celltype', 'sample', 'old', 'large', 'batch', 'sev.level', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'cell_type_original'
    uns: 'X_DR_sample', 'sample_embedding_params'
    obsm: 'X_DR_sample'
              sev.level  batch
sample                        
CoV-3-Wilk            4   Wilk
CoV-34-Aruna          4  Aruna
CoV-36-Aruna          4  Aruna
CoV-4-Wilk            4   Wilk
HD-1-Wilk             1   Wilk
HD-20-Aruna           1  Aruna
HD-21-Aruna           1  Aruna
HD-4-Wilk             1   Wilk
```

- `X` and `obsm['X_DR_sample']` are the embedding (8 samples × 7 PCs).
- `obs` takes each cell-level column's value from the first cell of the sample. That is right for sample-level columns such as `sev.level` and `batch`, but per-cell columns such as `total_counts` or `cell_type` then hold one arbitrary cell's value, not a per-sample summary.
- This object has no pseudotime. The pipeline adds `obs['pseudotime_sample']` only to its in-memory copy (`res["rna_results"]["pseudo_adata"]` in the [quick start](quickstart.md#3-the-same-run-from-python)); on disk it is in `CCA/pseudotime_sample.csv`.
- For multi-omics output, pass `modality_col="modality"` as well, to get one row per sample × modality.

## The sample embedding files

`sample_embedding/` keeps the embedding as a CSV and the four blocks it was built from:

```python
import numpy as np
import pandas as pd

emb = pd.read_csv(f"{OUT}/sample_embedding/sample_embedding.csv", index_col=0)
print(emb.round(2))
print("matches uns['X_DR_sample']:", np.allclose(emb.values, adata.uns["X_DR_sample"].values, atol=1e-5))

blocks = np.load(f"{OUT}/sample_embedding/sample_embedding_blocks.npz")
for key in blocks.files:
    print(f"{key:14s} {blocks[key].shape}")
```

```text
                PC1   PC2   PC3   PC4   PC5   PC6   PC7
sample                                                 
CoV-3-Wilk     4.41 -7.81 -5.21 -0.07  0.93  0.04  0.29
CoV-34-Aruna -12.39 -5.75  3.59  0.19 -0.15 -0.16 -0.01
CoV-36-Aruna  -1.62  2.35 -2.84  1.91 -3.84 -0.79 -0.53
CoV-4-Wilk    -5.90  5.42 -2.64 -5.03  1.00 -0.04  0.27
HD-1-Wilk      0.93  3.47  0.93  3.57  0.28  0.81  1.57
HD-20-Aruna    1.57  1.64  0.51  1.30  1.17  3.02 -1.16
HD-21-Aruna    3.34  2.09  1.42  1.75  2.25 -2.96 -0.61
HD-4-Wilk      9.65 -1.43  4.24 -3.64 -1.65  0.08  0.18
matches uns['X_DR_sample']: True
unit_ids       (8,)
A1_cell_types  (8,)
A1             (8, 8)
A2             (8, 120)
A3             (8, 299)
RMD            (8, 56)
```

These are the RNA results of the RNA + ATAC run (quick start step 6), so the numbers differ a little from steps 2–3: preprocessing and cell typing are not bit-for-bit reproducible between runs.

| Key in the `.npz` | Shape | Content |
| --- | --- | --- |
| `unit_ids` | (8,) | Sample IDs: the row order of every block and of the CSV. For multi-omics, `<sample>_<modality>`. |
| `A1_cell_types` | (8,) | The cell-type labels, i.e. the column order of `A1`. |
| `A1` | (8, 8) | Fraction of each sample's cells in each cell type. |
| `A2` | (8, 120) | Soft composition over k-means clusters of `Z_comp`, with `medium_K` = 120 clusters. |
| `A3` | (8, 299) | The same at `fine_K`, which asks for 300 but is capped at cells / 100 (299 here). |
| `RMD` | (8, 56) | For each cell type, the shift of the sample's mean in `Z_rmd` away from the mean of the other samples in the same batch (leave-one-out). PCA then reduces each cell type's shift to at most `rmd_dim` (8) dimensions and at most samples − 1 (7), so 8 cell types × 7 = 56. |

The blocks are saved before weighting. The embedding stage centres each block and scales it to a Frobenius norm of √(number of samples) × its weight in `block_weights` (next section). It then stacks the blocks side by side, runs PCA, and applies the sample-level Harmony step. The [sample embedding page](sample_embedding.md) covers the method in depth.

## Embedding parameters

`uns['sample_embedding_params']` records what the embedding stage actually used:

```python
for key, value in adata.uns["sample_embedding_params"].items():
    print(f"{key:18s} {value}")
```

```text
K_c                8
backend            cpu
batch_col          batch
batch_cols_multi   ['batch']
batch_method       harmony
block_weights      [6.11350963 1.5785014  1.         6.39      ]
cluster_emb_key    Z_comp
comp_emb_key       Z_comp
fine_K             299
medium_K           120
modality_col       
pca_components     10
rmd_emb_key        Z_rmd
rmd_weight         6.39
seed               42
use_clr            False
use_rmd            True
```

- `K_c`, `medium_K`, `fine_K`: the number of cell types and the k-means cluster counts actually used, after the caps.
- `block_weights`: the weights of A1, A2, A3 and RMD. The composition weights are √(`fine_K` / K): √(299/8) = 6.11 for A1, √(299/120) = 1.58 for A2 and 1 for A3.
- `rmd_weight`: the RMD weight α. The config value `"equal"` resolves to √(6.11² + 1.58² + 1²) = 6.39, which gives the RMD block as much weight as the three composition blocks together.
- `pca_components`: the number requested (`rna_sample_embedding_pca_components`). The embedding has at most samples − 1 PCs, so 7 here.
- `batch_col`, `batch_cols_multi`, `batch_method`: the sample-level batch correction that was requested. It ran here because there are 8 samples in 2 batches; with fewer than 8 samples, or a single batch, the PCA result is kept as is (the run log prints `[batch correction] 1 groups` in that case).
- `comp_emb_key`, `rmd_emb_key`: the cell embeddings that were used. `cluster_emb_key` repeats `comp_emb_key` under its old name.
- `backend`: `cpu` or `gpu`.

## Downstream result tables

All downstream tables are CSV files indexed by sample:

```python
dist = pd.read_csv(f"{OUT}/Sample_distance/cosine/distance_matrix_sample_DR.csv", index_col=0)
print(dist.iloc[:4, :4].round(2))

assoc = pd.read_csv(f"{OUT}/sample_association/variance_explained_sample.csv")
print(sorted(assoc["variable"].unique()))
print(assoc.query("variable == 'sev.level'")[["component", "r2", "perm_p", "fdr"]].round(3))
```

```text
              CoV-3-Wilk  CoV-34-Aruna  CoV-36-Aruna  CoV-4-Wilk
sample                                                          
CoV-3-Wilk          0.00          0.75          0.77        0.95
CoV-34-Aruna        0.75          0.00          0.64        0.48
CoV-36-Aruna        0.77          0.64          0.00        0.45
CoV-4-Wilk          0.95          0.48          0.45        0.00
['batch', 'cell_type', 'cell_type_original', 'celltype', 'large', 'old', 'pct_counts_mt', 'sev.level', 'total_counts', 'total_counts_mt']
   component     r2  perm_p    fdr
49       PC1  0.378    0.14  0.915
50       PC2  0.111    0.39  0.941
51       PC3  0.331    0.18  0.915
52       PC4  0.075    0.46  0.948
53       PC5  0.082    0.52  0.948
54       PC6  0.023    0.64  0.948
55       PC7  0.000    1.00  1.000
```

- The distance matrices are symmetric, samples × samples. The distances are rescaled so that the largest is 1.
- The association table has one row per PC × `obs` column: `r2` is the share of that PC's variance explained by the column, `perm_p` a permutation p-value (99 permutations, so never below 0.01), and `fdr` is adjusted over all rows. With 8 samples, `sev.level` explains 38% of PC1 and is still far from significant.
- With `rna_association_continuous_cols` and `rna_association_categorical_cols` left empty, every `obs` column is tested. That includes per-cell columns such as `total_counts` and `cell_type`, whose per-sample value is one cell's (see above). Ignore those rows, or list the sample-level columns you want in those two keys.

## The status file

Every run keeps a small JSON log, `sys_log/main_process_status.json`, recording which stages have finished for each pipeline. It is written when the run starts and again after each pipeline completes.

```python
import json

with open("sampledisco_demo_output/sys_log/main_process_status.json") as fh:
    status = json.load(fh)
print(sorted(status))
print(status["rna"])
```

```text
['atac', 'multiomics', 'rna', 'system_info']
{'preprocessing': True, 'cell_type_cluster': True, 'derive_sample_embedding': True, 'cca_based_cell_resolution_selection': False, 'sample_distance_calculation': True, 'trajectory_analysis': True, 'trajectory_dge': True, 'sample_cluster': True, 'proportion_test': True, 'cluster_dge': False, 'visualization': True, 'dimension_association_analysis': True}
```

- `cluster_dge` is `False` because RAISIN is off in the quick-start config. `cca_based_cell_resolution_selection` is a leftover that no current stage sets.
- `system_info` records the platform, the Python version, and whether a GPU was requested (`use_gpu`) and actually used (`gpu_available`).

The config key `initialization` decides what happens to this file at the start of a run:

- **`initialization: false`** (the demo config): the run reads the existing file, prints `Resuming from previous progress:` with its content, and merges its own results into it. A flag, once `True`, stays `True`, so the file records every stage that has ever finished in this `output_dir`, not only the last run's.
- **`initialization: true`** (the default when you call `wrapper()` without the key): the run resets every flag to `False`. It also deletes `<output_dir>/rna/result/`, `atac/result/` and `multiomics/result/` if they exist, but no current stage writes there, so in practice no earlier output is removed.

Neither value makes the pipeline skip a stage. Despite the "Resuming" message, which stages run is decided only by the stage flags in the config, as the next section shows. A rerun also never cleans up: files are overwritten in place, `trajectoryDEG/` gains a new timestamped set each time, and a folder from a stage you have since switched off stays where it was. To start clean, delete `output_dir` or point it to a new folder.

## Resuming a run

There is no "resume" switch. What a run does is set only by the stage flags in the config, and the early stages read their input from disk when you switch them off:

| To skip | Set | The run then |
| --- | --- | --- |
| QC and preprocessing | `rna_preprocessing: false` | reads `rna_adata_cell_path`, or, when that is empty, `<rna_output_dir>/preprocess/adata_preprocessed.h5ad`. |
| Cell typing | `rna_cell_type_cluster: false` | keeps the `cell_type` labels already in that file. |
| The sample embedding | `rna_derive_sample_embedding: false` | reuses `uns['X_DR_sample']` from that file. The run that computed it must have had `rna_save_cell_adata_after_embedding: true` (the template value). If the file has no embedding, the downstream stages are skipped with a message. |
| A downstream analysis | its own flag `false` | Each downstream flag is independent, except that trajectory DGE runs only together with the trajectory stage. |

The ATAC pipeline uses the same keys with the `atac_` prefix. Multi-omics resumes differently, through `multiomics_integration: false` and `multiomics_integrated_h5ad_path`.

!!! warning "Keys that do not resume anything"
    - `rna_adata_sample_path` and `rna_pseudo_adata_path`, and their `atac_` versions, are accepted by the config check but never read. The same goes for `multiomics_pseudobulk_h5ad_path`.
    - `initialization` and the status file never make the pipeline skip a stage (see [The status file](#the-status-file)).
    - `rna_sample_meta_path` and `rna_cell_meta_path` are read only by the preprocessing stage, so with `rna_preprocessing: false` they have no effect.

The three recipes below start from the quick-start `config.yaml` (RNA only) and build on each other. The first two take under a minute each on 8 cores (49 s each, about 3.5 GB of RAM), against almost 3 minutes for the full RNA run.

### Re-run only the downstream analysis

Copy the config, switch the first three stages off, and change something downstream, here the number of k-means clusters:

```bash
cp config.yaml config_resume.yaml
sed -i.bak \
    -e 's/^rna_preprocessing: true/rna_preprocessing: false/' \
    -e 's/^rna_cell_type_cluster: true/rna_cell_type_cluster: false/' \
    -e 's/^rna_derive_sample_embedding: true/rna_derive_sample_embedding: false/' \
    -e 's/^rna_cluster_number: 3/rna_cluster_number: 2/' \
    config_resume.yaml
sampledisco -m complex --config config_resume.yaml
cat sampledisco_demo_output/rna/sample_cluster/kmeans_clusters_sample.csv
```

```text
Resuming from previous progress:
{
  "rna": {
    "preprocessing": true,
    ...
}

============================================================
RUNNING RNA PIPELINE
============================================================
Starting RNA wrapper function...
RNA preprocessing pipeline completed successfully!
Starting sample distance calculation...
...
Total execution time: 33.93 seconds (0.57 minutes)
sample,cluster_sample_kmeans
CoV-3-Wilk,0
CoV-34-Aruna,1
CoV-36-Aruna,0
CoV-4-Wilk,1
HD-1-Wilk,0
HD-20-Aruna,0
HD-21-Aruna,0
HD-4-Wilk,0
```

The RNA wrapper prints no stage messages of its own this time: it only loads `preprocess/adata_preprocessed.h5ad`, with its cell types and embedding. Every downstream stage then reruns on that file. The clusters are now 0 and 1.

### Recompute the sample embedding, keep preprocessing and cell types

Switch the embedding stage back on and change one of its settings, here a fixed RMD weight α = 1 instead of `"equal"`:

```bash
sed -i.bak \
    -e 's/^rna_derive_sample_embedding: false/rna_derive_sample_embedding: true/' \
    -e 's/^rna_sample_embedding_rmd_weight: "equal"/rna_sample_embedding_rmd_weight: 1.0/' \
    config_resume.yaml
sampledisco -m complex --config config_resume.yaml
```

The log shows the new weights, with the RMD weight now 1.0 instead of 6.39, and the embedding written back into the h5ad:

```text
Starting sample embedding derivation (composition + RMD)...
[sample_embedding] weights=[6.114, 1.579, 1.0, 1.0] (n_blocks=4)
  [batch correction] 2 groups (group only)
[sample_embedding] updated ./sampledisco_demo_output/rna/preprocess/adata_preprocessed.h5ad (.uns['X_DR_sample'], in-place)
[sample_embedding] wrote ./sampledisco_demo_output/rna/sample_embedding/sample_embedding.csv
[sample_embedding] wrote ./sampledisco_demo_output/rna/sample_embedding/sample_embedding_blocks.npz
```

Every downstream stage then reruns on the new embedding. Because the h5ad now holds this embedding, a later run with `rna_derive_sample_embedding: false` reuses it.

### Add a stage later: RAISIN

RAISIN is off in the demo config because of its memory use. To add it to a finished run, switch the embedding stage off again and RAISIN on:

```bash
sed -i.bak \
    -e 's/^rna_derive_sample_embedding: true/rna_derive_sample_embedding: false/' \
    -e 's/^rna_cluster_dge: false/rna_cluster_dge: true/' \
    config_resume.yaml
sampledisco -m complex --config config_resume.yaml
ls sampledisco_demo_output/rna/raisin_results/*
```

```text
sampledisco_demo_output/rna/raisin_results/1_vs_4:
raisin_results.csv
volcano_labeled.png
volcano.png

sampledisco_demo_output/rna/raisin_results/summary_plots:
all_results_combined.csv
cluster_gene_zscore.png
pseudobulk_heatmap.png
raisin_summary.csv
raisin_summary.txt
summary_dotplot.png
```

- On 8 CPU cores this took 2 min 58 s and peaked at 15.7 GB of RAM, against about 3.5 GB for the other two recipes. The ComBat batch correction inside RAISIN works on a dense genes × cells matrix. Make sure the machine has the memory before you switch it on for a larger dataset.
- `1_vs_4` compares the two `sev.level` groups, not the k-means clusters: `rna_cluster_differential_gene_group_col: "sev.level"` in the demo config takes precedence over the clusters, for RAISIN and for the proportion test. Set that key to empty (`~`) to compare clusters instead. With 8 samples some clusters hold a single sample, and from 0.6.0 RAISIN skips every comparison involving such a group with a warning and a `FAILED` marker file. A one-sample group has no residual degrees of freedom; 0.5.0 set its variance to 0, which produced spurious "significant" genes.
- RAISIN starts one thread per usable CPU: the CPUs the process may run on, capped at `SLURM_CPUS_PER_TASK` under Slurm. The log says `Using 8 CPU cores` on this 8-CPU job. Before 0.6.0 it counted both hardware threads of every core (16 here).
- The other downstream stages rerun as well. Switch their flags off if you only want RAISIN.
