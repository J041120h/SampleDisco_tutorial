# Outputs and resuming runs

A config run writes everything under `output_dir`: one folder per pipeline (`rna/`, `atac/`, `multiomics/`) and a `sys_log/` folder. The snippets below run from the folder of the RNA-only [Quick start](quickstart.md) run, whose `output_dir` is `./sampledisco_demo_output`.

## The output tree

```text
sampledisco_demo_output/
├── sys_log/main_process_status.json      which stages have run
├── rna/
│   ├── preprocess/
│   │   ├── adata_preprocessed.h5ad       the cell-level AnnData (200 MB on the demo)
│   │   ├── cell_type.csv
│   │   └── umap_cell_type.png
│   ├── sample_embedding/
│   │   ├── sample_embedding.csv          samples × PCs
│   │   └── sample_embedding_blocks.npz   the four blocks, before weighting
│   ├── Sample_distance/<method>/
│   ├── CCA/                              TSCAN/ for an unsupervised trajectory
│   ├── trajectoryDEG/
│   ├── sample_cluster/                   k-means clusters, with proportion_test/ inside
│   ├── raisin_results/                   only with rna_cluster_dge: true
│   ├── sample_association/
│   └── visualization/
├── atac/                                 with run_atac_pipeline: true; the same layout as rna/
└── multiomics/                           with run_multiomics_pipeline: true; adata_sample.h5ad in preprocess/
```

A rerun overwrites files in place and never cleans up: the folder of a stage you have since switched off stays. To start clean, delete `output_dir` or use a new one.

## Stage by stage

| Stage (flag) | Folder | Files |
|---|---|---|
| Preprocessing and cell typing (`*_preprocessing`, `*_cell_type_cluster`) | `preprocess/` | `adata_preprocessed.h5ad` (see below), `cell_type.csv`, `umap_cell_type.png` |
| Sample embedding (`*_derive_sample_embedding`) | `sample_embedding/` | `sample_embedding.csv`, `sample_embedding_blocks.npz` ([what they hold](sample_embedding.md#outputs)). Also stored in the h5ad, unless `*_save_cell_adata_after_embedding: false`. |
| Sample distance (`*_sample_distance_calculation`) | `Sample_distance/<method>/` | `distance_matrix_sample_DR.csv`, a heatmap PDF, `distance_check_results_sample_DR_<method>.txt` (mean distance within and between groups) |
| Trajectory (`*_trajectory_analysis`) | `CCA/` or `TSCAN/` | CCA: `pseudotime_sample.csv` and two PDFs. TSCAN: `X_DR_sample_pseudotime.csv` and two cluster plots. |
| Trajectory DGE (`*_trajectory_dge`) | `trajectoryDEG/` | `gam_all_genes_<timestamp>.tsv`, `gam_significant_<timestamp>.tsv`, `gam_summary_<timestamp>.txt`, `visualizations/`. Each run adds a new timestamped set. |
| Sample clustering (`*_sample_cluster`) | `sample_cluster/` | `kmeans_clusters_sample.csv`, `kmeans_sample_embedding.png` |
| Proportion test (`*_proportion_test`) | `sample_cluster/proportion_test/` | `proportion_test_<a>_vs_<b>.csv` per pair of groups, `proportion_test_significant_summary.txt`, heatmaps |
| RAISIN (`*_cluster_dge`) | `raisin_results/` | `<a>_vs_<b>/raisin_results.csv` with volcano plots, and `summary_plots/` |
| Visualization (`*_visualize_data`) | `visualization/` | `cell_type_dendrogram.png`, and `proportion_embedding_<col>.png` and `expression_embedding_<col>.png` per grouping column |
| Dimension association (`*_dimension_association_analysis`) | `sample_association/` | `variance_explained_sample.csv`, `figures/` |

Multi-omics keeps its cell-level object, with the cells of both modalities, in `preprocess/adata_sample.h5ad`, and draws `visualization/all_samples_*_by_<col>.png`. A run that trains scGLUE also writes `integration/glue/` and `preprocess/adata_{rna,atac}_preprocessed.h5ad`. What each analysis computes is on the [Downstream analysis](downstream/index.md) page.

## Load the results back

The cell-level file holds everything up to the sample embedding. `build_sample_adata` turns it into the one-row-per-sample AnnData that every downstream function takes:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

OUT = "sampledisco_demo_output/rna"
adata = ad.read_h5ad(f"{OUT}/preprocess/adata_preprocessed.h5ad")
pseudo = build_sample_adata(adata, sample_col="sample")
print(pseudo)
```

```text
[build_sample_adata] dropped 10 obs column(s) that vary within a sample (per-cell values): ['celltype', 'old', 'large', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'cell_type_original']
AnnData object with n_obs × n_vars = 8 × 7
    obs: 'sample', 'batch', 'sev.level'
    uns: 'X_DR_sample', 'sample_embedding_params'
    obsm: 'X_DR_sample'
```

| In `adata` | Content |
|---|---|
| `X`, `layers['counts']` | log-normalized expression of the genes that passed QC, and the raw counts |
| `obs['cell_type']` | the labels from cell typing (`cell_type_original`: the Leiden clusters before merging) |
| `obsm['Z_comp']`, `obsm['Z_rmd']` | the cell embeddings with sample differences removed and kept ([legacy names](../legacy_keys.md)) |
| `uns['X_DR_sample']` | the sample embedding, a DataFrame of samples × PCs |
| `uns['sample_embedding_params']` | the settings it was computed with ([details](sample_embedding.md#outputs)) |

- `pseudo` has the embedding in `X` and `obsm['X_DR_sample']`, and only the sample-level `obs` columns: those constant within each sample. For multi-omics output, pass `modality_col="modality"` as well, to get one row per sample × modality.
- The downstream results, including the pseudotime, are CSV files indexed by sample, e.g. `pd.read_csv(f"{OUT}/CCA/pseudotime_sample.csv")`.

## The status file

`sys_log/main_process_status.json` records which stages have finished for each pipeline:

```python
import json

with open("sampledisco_demo_output/sys_log/main_process_status.json") as fh:
    status = json.load(fh)
print(status["rna"])
```

```text
{'preprocessing': True, 'cell_type_cluster': True, 'derive_sample_embedding': True, 'cca_based_cell_resolution_selection': False, 'sample_distance_calculation': True, 'trajectory_analysis': True, 'trajectory_dge': True, 'sample_cluster': True, 'proportion_test': True, 'cluster_dge': False, 'visualization': True, 'dimension_association_analysis': True}
```

- `status["system_info"]` records the platform, whether a GPU was requested (`use_gpu`) and whether it was used (`gpu_available`).
- `initialization: false` (the template) reads the file, prints `Resuming from previous progress:` and merges the new run into it: a flag, once `True`, stays `True`. `initialization: true` (the `wrapper()` default) resets every flag to `False`.
- **Neither value skips a stage.** Despite the "Resuming" message, only the stage flags in the config decide what runs.

## Resuming a run

There is no "resume" switch. You switch finished stages **off**, and the early stages then read their results from disk:

| To skip | Set | The run then |
|---|---|---|
| QC and preprocessing | `rna_preprocessing: false` | reads `rna_adata_cell_path`, or by default `{rna output}/preprocess/adata_preprocessed.h5ad` |
| Cell typing | `rna_cell_type_cluster: false` | keeps the labels already in that file (`rna_celltype_col`, normally `cell_type`) |
| The sample embedding | `rna_derive_sample_embedding: false` | reuses `uns['X_DR_sample']` from that file; without it every downstream stage is skipped |
| A downstream analysis | its own flag `false` | Downstream stages are independent and rerun from the embedding, except that trajectory DGE runs only with the trajectory stage. |

- ATAC uses the same keys with `atac_`. Multi-omics resumes through `multiomics_integration: false` ([Configuration guide](configuration.md#multi-omics)).
- Keeping `*_cell_type_cluster: true` re-types the cells and overwrites `adata_preprocessed.h5ad`. `*_sample_meta_path` and `*_cell_meta_path` are read only by preprocessing.
- **Keys that do not resume anything:** `rna_adata_sample_path`, `rna_pseudo_adata_path`, their `atac_` versions and `multiomics_pseudobulk_h5ad_path` are required by the CLI but never read. Leave them `null`.

### Re-run only the downstream analysis

Switch the first three stages off and change something downstream, here the number of k-means clusters:

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

The run loads `adata_preprocessed.h5ad` with its cell types and embedding and reruns every downstream stage. On 8 CPU cores it took under a minute. `kmeans_clusters_sample.csv` now holds two clusters, which here separate the healthy (`HD-*`) from the COVID-19 (`CoV-*`) samples.

### Add a stage later: RAISIN

RAISIN is off in the template because of its memory use. To add it to a finished run, switch it on in the same config:

```bash
sed -i.bak 's/^rna_cluster_dge: false/rna_cluster_dge: true/' config_resume.yaml
sampledisco -m complex --config config_resume.yaml
ls sampledisco_demo_output/rna/raisin_results/
```

```text
1_vs_4
summary_plots
```

- On 8 CPU cores this peaked at about 16 GB of RAM, against 3.5 GB without RAISIN. The other downstream stages reran too (about 5 minutes in all, mostly RAISIN); switch their flags off if you only want RAISIN.
- `1_vs_4` compares the two `sev.level` groups, because the template sets `rna_cluster_differential_gene_group_col: "sev.level"`. Set that key to `null` to compare the k-means clusters instead.
