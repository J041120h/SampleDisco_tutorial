# Multi-omics from a config

This page runs the RNA + ATAC branch with one command, `sampledisco --config`, on the COVID-19 demo data. It does what the [step-by-step multi-omics page](multiomics.md) does, then the downstream analyses:

1. Start from the pre-integrated scGLUE object (no training).
2. Cell typing, sample embedding and embedding plots.
3. Sample distances, CCA trajectory, sample clustering, proportion test and dimension association.

The last section trains scGLUE from the raw counts instead.

You need a working installation (see [Installation](../installation.md)). The pre-integrated route needs no scGLUE.

!!! note "Run everything from one folder"
    The config uses relative paths. Run every command on this page from the same folder. It can be the folder of the [step-by-step page](multiomics.md): that page writes to `sampledisco_demo_output/multiomics_api/`, this one to `sampledisco_demo_output/multiomics/`.

## 1. Set up a run folder

```bash
mkdir -p data
curl -fL -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
```

`test_multiomics_integrated.h5ad` (54 MB) holds the scGLUE embedding of 29,989 RNA and 29,021 ATAC cells from 8 + 8 donors. See [what it contains](multiomics.md#inputs). Keep the `-f`: without it, curl saves an HTTP error page under the `.h5ad` name and still reports success.

Write the demo config:

```bash
sampledisco --init-config config.yaml
```

`--init-config` writes the packaged demo config, the same file as [`config_demo.yaml`](../assets/config_demo.yaml). It refuses to overwrite an existing file.

## 2. Edit the config

The template runs RNA and ATAC and leaves multi-omics off. Switch the single-modality pipelines off and multi-omics on, and color the embedding plots by `sev.level`. The template already points `multiomics_integrated_h5ad_path` at `./data/test_multiomics_integrated.h5ad`.

```bash
sed -i.bak \
    -e 's/^run_rna_pipeline: true/run_rna_pipeline: false/' \
    -e 's/^run_atac_pipeline: true/run_atac_pipeline: false/' \
    -e 's/^run_multiomics_pipeline: false.*/run_multiomics_pipeline: true/' \
    -e 's/^multiomics_visualization_grouping_column:.*/multiomics_visualization_grouping_column: ["sev.level"]/' \
    config.yaml
grep -E '^(output_dir|use_gpu|run_[a-z]+_pipeline|multiomics_(integration|integrated_h5ad_path|use_gpu|visualization_grouping_column|trajectory_dge|cluster_dge)):' config.yaml
```

```text
output_dir: "./sampledisco_demo_output"
run_rna_pipeline: false
run_atac_pipeline: false
run_multiomics_pipeline: true
use_gpu: false   # set true to use the RAPIDS GPU paths if installed; auto-falls back to CPU otherwise
multiomics_integration: false   # false -> SKIP scGLUE training, load the pre-integrated file below. Set true to train from scratch (needs bedtools) + set the resume path to null.
multiomics_trajectory_dge: false   # 16 (sample, modality) units: GAM cannot fit on small test
multiomics_cluster_dge: false      # disabled: same RAISIN OOM risk as ATAC (large peak matrix)
multiomics_integrated_h5ad_path: "./data/test_multiomics_integrated.h5ad"   # scGLUE-integrated object from Zenodo; used when multiomics_integration: false (skips training)
multiomics_use_gpu: true
multiomics_visualization_grouping_column: ["sev.level"]
```

`sed -i.bak` works with both GNU and macOS `sed` and keeps the original as `config.yaml.bak`. You can just as well change the lines in a text editor. The keys that matter here:

| Key | Value | What it does |
| --- | --- | --- |
| `run_multiomics_pipeline` | `true` | runs the multi-omics branch |
| `multiomics_integration` | `false` (template) | skips scGLUE and starts from the integrated object instead |
| `multiomics_integrated_h5ad_path` | `./data/test_multiomics_integrated.h5ad` (template) | the integrated object; read only, never written |
| `multiomics_visualization_grouping_column` | `["sev.level"]` | colors the embedding plots; with the template's empty value they are drawn uncolored |
| `use_gpu`, `multiomics_use_gpu` | `false`, `true` (template) | the run stays on CPU; see [GPU](#gpu) |
| `multiomics_trajectory_dge`, `multiomics_cluster_dge` | `false` (template) | cannot run on this route and are skipped if switched on; see [DGE and RAISIN](#dge-and-raisin-need-the-per-modality-rna-file) |

!!! note "Where the results go"
    The input file is not modified. Every update goes to `<multiomics_output_dir>/preprocess/adata_sample.h5ad`, here `sampledisco_demo_output/multiomics/preprocess/adata_sample.h5ad`: the run starts it as a full copy of the input (57 MB) and adds the cell types and the sample embedding. A from-scratch run (`multiomics_integration: true`) writes its object to the same file.

    The pipeline reads that output file first and falls back to `multiomics_integrated_h5ad_path` only when it does not exist. A second run in the same folder therefore starts from the first run's results; delete the output file to start again from the input.

!!! warning "sampledisco 0.5.0 wrote the results into the input file"
    In 0.5.0 the run writes the new cell types (cells of dropped clusters removed) and the sample embedding into `multiomics_integrated_h5ad_path` itself, and the output file above is never created. With 0.5.0, point the key at a copy in that location before running, so that the download stays untouched and the next steps find the results:

    ```bash
    mkdir -p sampledisco_demo_output/multiomics/preprocess
    cp data/test_multiomics_integrated.h5ad sampledisco_demo_output/multiomics/preprocess/adata_sample.h5ad
    sed -i.bak 's|^multiomics_integrated_h5ad_path:.*|multiomics_integrated_h5ad_path: "./sampledisco_demo_output/multiomics/preprocess/adata_sample.h5ad"|' config.yaml
    ```

    0.6.0 accepts this setup too and gives the same results.

## 3. Run the pipeline

```bash
sampledisco --config config.yaml
```

- `-m complex` is the default mode, so `sampledisco -m complex --config config.yaml` does the same.
- The CLI checks the config against `wrapper()` before it starts: every key must be present and unknown keys are rejected. Change values; don't delete keys.
- On 8 CPU cores the run took 2 to 5.5 minutes, depending on the node, and peaked at 2.8 GB of RAM. scGLUE is not needed.

The log walks through these stages:

1. `GPU mode: Disabled`, then the warning about the per-modality files (see [below](#dge-and-raisin-need-the-per-modality-rna-file)).
2. `Loaded embedding union from: ./data/test_multiomics_integrated.h5ad`, then `Updates are written to ./sampledisco_demo_output/multiomics/preprocess/adata_sample.h5ad; the input file is not modified.` After the run, `md5sum data/test_multiomics_integrated.h5ad` still prints the value in [Demo data](demo_data.md).
3. `Step 2b: obsm['Z_clust'] already present … — no Harmony pass needed.` The file already has a sample-removed embedding (under its old name), so `harmonize_xglue` is skipped. The message's "(end-to-end from 2-run scGLUE)" only names one possible source; here it was the older run.
4. `Step 2c: Cell type assignment (use_rep=Z_clust)`: 16 clusters, none dropped. The pipeline reads the old key names itself, with a `FutureWarning` for each.
5. `Step 3: Sample embedding`: `comp_emb=Z_clust, rmd_emb=Z_cmd`, 16 units × 10 PCs.
6. The downstream steps: sample distances (cosine and correlation), CCA trajectory on `sev.level`, k-means sample clustering (3 clusters), the proportion test between the `sev.level` groups, the embedding plots and the dimension association. The run ends with `Multiomics pipeline completed successfully!`.

## 4. What it writes

Everything goes to `sampledisco_demo_output/multiomics/`; `data/` is not touched.

| Path | Content |
| --- | --- |
| `preprocess/adata_sample.h5ad` | a copy of the input with the new `obs['cell_type']`, `obs['label_transfer_confidence']`, a new `X_umap`, and `uns['X_DR_sample']` + `uns['sample_embedding_params']` |
| `preprocess/cell_type.csv`, `preprocess/visualization/` | cell types per cell; the cell-type UMAPs, `label_transfer_confidence.png`, `celltype_modality_heatmap.png` and its CSV (as on the [step-by-step page](multiomics.md#2-joint-cell-typing)) |
| `sample_embedding/sample_embedding.csv`, `sample_embedding_blocks.npz` | the 16 × 10 sample embedding and its four blocks |
| `visualization/all_samples_{expression,proportion,combined}_by_sev.level.png` | the embedding colored by `sev.level` |
| `Sample_distance/{cosine,correlation}/` | unit × unit distance matrix, heatmap PDF and the in-group / between-group check |
| `CCA/` | `pseudotime_sample.csv` and the `pca_5d_cca_sample*.pdf` plots |
| `sample_cluster/` | `kmeans_clusters_sample.csv`, `kmeans_sample_embedding.png` and `proportion_test/` |
| `sample_association/` | `variance_explained_sample.csv` and two figures |
| `../sys_log/main_process_status.json` | which stages ran |

See [Understanding the outputs](outputs.md) for how to read each file, and the [Downstream analysis tutorials](downstream/index.md) for the analyses themselves.

What the 16 units show here:

- **Distances:** units of the same `sev.level` are closer than units of different levels (cosine: mean 0.48 within a group, 0.64 between groups; score 1.33).
- **CCA:** the best pair of PCs for `sev.level` is PC1 + PC5 (score 0.853).
- **Association:** `sev.level` explains 46% of the variance of PC1 (permutation p = 0.01, the smallest possible with the demo's 99 permutations); `modality` explains 19% of PC1 (p = 0.11).
- **Proportion test:** no cell type differs between the groups at FDR < 0.01.

!!! note "Small-sample caveats"
    - With 8 donors per modality these numbers are a sanity check, not a finding. The CCA score is the best of all PC pairs, so it is optimistic; `multiomics_cca_pvalue` stays `false` on the demo.
    - The association step also tests per-cell columns such as `n_genes`, `total_counts`, `cell_type` and `label_transfer_confidence`, because it takes every `obs` column, using the value of the first cell of each unit. Read only the rows of real sample-level variables (`sev.level`, `modality`), or set `multiomics_association_categorical_cols` and `multiomics_association_continuous_cols` yourself.
    - For a cohort large enough to test, see the [Sound Life 32-donor tutorial](soundlife_rna.md).

## 5. Load the results in Python

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

OUT = "sampledisco_demo_output/multiomics"
adata = ad.read_h5ad(f"{OUT}/preprocess/adata_sample.h5ad")
params = adata.uns["sample_embedding_params"]
print(adata.n_obs, "cells;", params["comp_emb_key"], "/", params["rmd_emb_key"], "on", params["backend"])

sample_adata = build_sample_adata(adata, sample_col="sample", modality_col="modality")
print(sample_adata)
```

```text
59010 cells; Z_clust / Z_cmd on cpu
AnnData object with n_obs × n_vars = 16 × 10
    obs: 'sample', 'sev.level', 'balancing_weight', 'modality', 'original_barcode', 'original_sample', 'n_genes_by_counts', 'total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'pct_counts_mt', 'n_genes', 'cell_type', 'label_transfer_confidence'
    uns: 'X_DR_sample', 'sample_embedding_params'
    obsm: 'X_DR_sample'
```

`sample_embedding_params` records which embeddings were used: the legacy `Z_clust` / `Z_cmd`, which hold the same values as the `Z_comp` / `Z_rmd` of the step-by-step page. `sample_adata` is the one-row-per-unit object that every downstream function accepts.

??? note "Same result as the Python route"
    If you also ran the [step-by-step page](multiomics.md) in this folder, compare the two sample embeddings:

    ```python
    import pandas as pd

    api = pd.read_csv("sampledisco_demo_output/multiomics_api/sample_embedding/sample_embedding.csv", index_col=0)
    cfg = pd.read_csv(f"{OUT}/sample_embedding/sample_embedding.csv", index_col=0)
    print("max |difference|:", (api - cfg.loc[api.index]).abs().to_numpy().max())
    ```

    ```text
    max |difference|: 0.0
    ```

    The two routes give the same cell types and the same sample embedding. Run on different machines, the numbers can differ in the fifth decimal (floating-point arithmetic); the cell types stay identical.

## DGE and RAISIN need the per-modality RNA file

The integrated object holds embeddings only: its `.X` has no expression values. Trajectory DGE and RAISIN therefore read the RNA cells from `<multiomics_output_dir>/preprocess/adata_rna_preprocessed.h5ad` instead. They use the RNA units of the sample embedding (`CoV-3-Wilk_RNA` becomes `CoV-3-Wilk`).

That file is written only when the run integrates the data itself (`multiomics_integration: true`, step `multiomics_run_glue_preprocess_per_modality`). With `multiomics_integration: false` it does not exist, so the two steps cannot run on this route, and the template leaves them off. The log says so at the start:

```text
WARNING: integration=False — per-modality preprocessed h5ads (./sampledisco_demo_output/multiomics/preprocess/adata_rna_preprocessed.h5ad, ./sampledisco_demo_output/multiomics/preprocess/adata_atac_preprocessed.h5ad) will NOT be (re)generated. propagate_cell_type and downstream DGE/RAISIN steps that depend on these files will be silently skipped unless they already exist on disk at those paths.
```

If you switch `multiomics_trajectory_dge` or `multiomics_cluster_dge` on anyway, the pipeline skips both steps with one more warning, runs every other downstream step, and finishes normally. We tested this with both keys on (and `multiomics_spline_order: 2`):

```text
WARNING: skipping multi-omics trajectory DGE and cluster DGE (RAISIN): the per-modality RNA expression file ./sampledisco_demo_output/multiomics/preprocess/adata_rna_preprocessed.h5ad was not found (it is written when multiomics_integration is true). The other downstream steps still run.
```

sampledisco 0.5.0 did not skip them: the whole run stopped with `RuntimeError: No features survived pseudobulk aggregation`.

To get DGE and RAISIN for multi-omics, [train scGLUE from the config](#train-scglue-from-the-config), which writes the per-modality files.

## GPU

`multiomics_use_gpu: true` alone does nothing. The pipeline uses the GPU only when **both** `use_gpu` and `multiomics_use_gpu` are `true`, on Linux. The template has `use_gpu: false`, so this run stays on CPU even though `multiomics_use_gpu` is `true`; the log says `GPU mode: Disabled`. With both on, `harmonize_xglue`, the cell typing (RAPIDS) and the sample embedding (cupy/cuML) run on the GPU, each falling back to the CPU when its libraries are missing. See [GPU acceleration](gpu.md#8-multi-omics-on-the-gpu) for a GPU run of this route.

scGLUE training is different: it uses a visible CUDA GPU on its own, whatever these two keys say.

## Train scGLUE from the config

To build the integration from the raw counts, download them and turn integration on. scGLUE, `bedtools` and the gene annotation must be available; see the [requirements](multiomics.md#requirements).

```bash
curl -fL -o data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
```

```bash
sed -i.bak \
    -e 's/^multiomics_integration: false.*/multiomics_integration: true/' \
    -e 's|^multiomics_integrated_h5ad_path:.*|multiomics_integrated_h5ad_path:|' \
    -e 's|^multiomics_output_dir:.*|multiomics_output_dir: "./sampledisco_demo_output/multiomics_glue_config"|' \
    -e 's/^multiomics_cluster_dge: false.*/multiomics_cluster_dge: true/' \
    config.yaml
sampledisco --config config.yaml
```

- The run does what [step 5 of the Python page](multiomics.md#5-train-scglue-from-scratch-optional) does: `multiomics_preparation` (`Step 1: Running GLUE integration...`), then `harmonize_xglue` (`Step 2b: Harmony post-pass on obsm['X_glue'] → obsm['Z_comp']`), cell typing, `propagate_cell_type`, the sample embedding and the downstream analyses.
- Emptying `multiomics_integrated_h5ad_path` keeps the pre-integrated file out of this run, and the new `multiomics_output_dir` keeps its results apart from section 3's.
- `multiomics_cluster_dge: true` now works: RAISIN reads the RNA cells from `preprocess/adata_rna_preprocessed.h5ad` and compares the `sev.level` groups (`multiomics_cluster_differential_gene_group_col: "sev.level"`). In our run it tested 9,856 genes and found none significant at FDR < 0.05, which is expected with 4 + 4 donors.
- On an NVIDIA L40S GPU node with 8 CPU cores the whole run took 30 minutes and peaked at 14.6 GB of RAM: 19 minutes for the scGLUE step (including the first-time annotation download; see the [requirements](multiomics.md#requirements)), 8 minutes for RAISIN (1 minute in another of our environments), and about 3 minutes for the rest.

**Writes** → everything listed in [section 4](#4-what-it-writes), under `sampledisco_demo_output/multiomics_glue_config/`, plus `integration/glue/` (the scGLUE model and inputs), `preprocess/adata_{rna,atac}_preprocessed.h5ad`, `visualization/scglue_umap_modality.png` and `raisin_results/` (`1_vs_4/raisin_results.csv`, volcano plots and `summary_plots/`).

The scGLUE keys, with the demo values:

| Key | Demo | Meaning |
| --- | --- | --- |
| `multiomics_rna_file`, `multiomics_atac_file` | `./data/test_RNA.h5ad`, `./data/test_ATAC.h5ad` | raw counts; `.obs` must have the sample column |
| `multiomics_ensembl_release`, `multiomics_species` | `98`, `homo_sapiens` | gene annotation for the guidance graph, downloaded by pyensembl on first use |
| `multiomics_n_top_genes`, `multiomics_n_top_peaks` | `2000`, `50000` | highly variable genes and peaks given to scGLUE |
| `multiomics_n_pca_comps`, `multiomics_n_lsi_comps` | `50`, `50` | RNA PCA and ATAC LSI dimensions |
| `multiomics_batch_col` | empty | a technical batch column. scGLUE removes it only when both files have it (the demo ATAC file has none); it is also added to the `harmonize_xglue` pass and the sample-level Harmony |
| `multiomics_glue_data_batch_size` | `128` | cells per minibatch (the function default is 1024) |
| `multiomics_glue_max_epochs` | `50` | training epochs cap; empty lets scGLUE choose |
| `multiomics_glue_dataloader_num_workers` and the other `multiomics_glue_*` workers | `0` | data-loading workers; the template keeps them at 0 (its comment cites a data-loader race on CPU) |
| `multiomics_harmonize_xglue_max_iter` | `30` | Harmony iterations for `Z_comp` |
| `multiomics_run_glue_twice_for_sample_removal` | `false` | `true` trains scGLUE a second time to get `Z_comp` instead of using Harmony |
| `multiomics_consistency_threshold` | `0.05` | only prints whether scGLUE's integration-consistency score passes; it changes nothing |

---

Continue with the [Downstream analysis tutorials](downstream/index.md) for what each downstream output means, and with the [Configuration guide](configuration.md) for the full key reference.
