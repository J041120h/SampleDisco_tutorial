# Configuration guide

A SampleDisco run is one YAML file. The CLI reads it into a flat dictionary and calls [`wrapper(**config)`](../api/pipeline/wrapper.md): every key is a parameter of `wrapper()`. This page explains how the keys fit together: which stage each one drives, which stages depend on others, how to resume a run, and which keys look meaningful but do nothing. Every key, with its default and the function that consumes it, is listed in the generated [Config key reference](config_reference.md).

If you have not run SampleDisco yet, do the [Quick start](quickstart.md) first. To point a config at your own data, read the box below and then [Using your own data](own_data.md).

!!! warning "Change these keys before you run on your own data"
    The template and the `wrapper()` defaults are wired to the COVID-19 demo, whose phenotype column is `sev.level`. Six keys per modality name that column, and **`null` does not switch them off: it falls back to `sev.level`.**

    | Key (each with an `rna_`, `atac_` and `multiomics_` version) | Demo value | What happens if the column is missing |
    |---|---|---|
    | `*_grouping_columns` | `["sev.level"]` | The sample-distance step skips its group-separation score (`Skipping group distance check …`), then the visualization step raises `ValueError: Missing required columns in pseudobulk AnnData: ['sev.level']` and, with `continue_on_error: false`, the run stops. |
    | `*_trajectory_col` | `"sev.level"` | With `*_trajectory_supervised: true`, the trajectory step raises `ValueError: Trajectory column 'sev.level' not found in pseudo_adata.obs.` With TSCAN it only orients the trajectory-DGE pseudotime (a warning if missing). |
    | `*_trajectory_visualization_label` | `["sev.level"]` | TSCAN only: its group plot is skipped with a warning. A numeric column here also picks the start of the trajectory when `*_tscan_origin` is `null`. |
    | `*_cluster_differential_gene_group_col` | `"sev.level"` | The proportion test catches the `KeyError`, prints it and produces nothing; RAISIN falls back to the k-means clusters. |
    | `*_association_categorical_cols` / `*_association_continuous_cols` | `null` (every sample-level column) | Nothing breaks; set them to the covariates you care about. |
    | `*_autotune_grouping_col` | `null` | Only read with `*_autotune_enable: true`. |

    Point them at a sample-level column of your data (one value per sample). `*_trajectory_col` must be **numeric** for the supervised (CCA) trajectory: text values become NaN and those samples are dropped. Also check `*_min_cells`, `*_min_genes` and the other QC thresholds, which are tuned for the demo; `sampledisco --init-config` prints the same reminder.

## 1. Get a template

Start from the template that ships inside the package. It is complete, every key is present, and it runs on the demo data as written.

```bash
sampledisco --init-config my_config.yaml
```

```text
Wrote the demo config to 'my_config.yaml'.
Edit the data paths / options, then run:
  sampledisco -m complex --config my_config.yaml
WARNING: 'my_config.yaml' is the config used for the Zenodo demo dataset — its QC thresholds and trajectory/grouping columns (e.g. 'sev.level', '*_supervised: true') are tuned for that demo data, not yours. Review and adjust every field before running on your own data.
```

- Without a path, `--init-config` writes `./config.yaml`. It refuses to overwrite an existing file.
- The file is identical to [`config_demo.yaml`](../assets/config_demo.yaml). It reads `./data/test_RNA.h5ad` and `./data/test_ATAC.h5ad` (see [Demo data](demo_data.md)) and writes to `./sampledisco_demo_output/`.
- Edit values; don't delete keys. The CLI rejects a config with a missing key.

## 2. Run it

=== "CLI"

    ```bash
    sampledisco -m complex --config my_config.yaml
    ```

    - `-m complex` is the only mode and the default; `python -m sampledisco.cli -m complex --config my_config.yaml` is equivalent (`python -m sampledisco` does not work).
    - Before anything runs, the CLI checks the file against `wrapper()`: an unknown key raises `ValueError: Unexpected parameter in config: '<key>'`, a missing one `ValueError: Missing required parameter in config: '<key>'`. Only the three `*_save_cell_adata_after_embedding` keys may be left out. `*_sample_embedding_rmd_weight` is also checked (`"equal"` or a number ≥ 0).
    - The exit status is 1 if any pipeline failed, including with `continue_on_error: true`.

=== "Python"

    ```python
    import yaml
    from sampledisco.wrapper.wrapper import wrapper

    with open("my_config.yaml") as fh:
        config = yaml.safe_load(fh)

    results = wrapper(**config)
    print(sorted(results))
    ```

    ```text
    ['atac_results', 'rna_results', 'status_flags', 'system_info']
    ```

    - Import from `sampledisco.wrapper.wrapper`: it works on 0.5.0 and 0.6.0. `from sampledisco import wrapper` also works on 0.6.0; in 0.5.0 it gives a module, and calling it raises `TypeError: 'module' object is not callable`.
    - `wrapper()` does **not** validate the key set. A misspelled key raises `TypeError: wrapper() got an unexpected keyword argument`, but a missing key silently takes its default, which is often not the demo value (see [section 10](#10-python-with-a-partial-config)).
    - `results["rna_results"]["pseudo_adata"]` is the sample-level `AnnData` (one row per sample, embedding in `obsm['X_DR_sample']`), and `results["rna_results"]["adata"]` the cell-level one. Same for `atac_results` and `multiomics_results`.

On the demo template (RNA and ATAC, 8 CPU cores) a run took between 3 min and 8 min 40 s depending on the node, with a peak of 7.1 GB RAM. Running it a second time with `initialization: false` recomputed every stage again: nothing is skipped (see [section 9](#9-resuming-a-run)).

## 3. How keys map to functions

The config is flat. The prefix says which pipeline reads a key; the rest of the name is usually the parameter of the function it ends up in.

- `rna_*` keys go to the RNA pipeline, `atac_*` keys to the ATAC pipeline and `multiomics_*` keys to the multi-omics pipeline.
- Keys without a prefix are global (section 5).

| Config key | Becomes | In function |
|---|---|---|
| `rna_min_genes: 200` | `min_genes=200` | `preprocess` (RNA) |
| `atac_leiden_cluster_resolution: 0.8` | `leiden_cluster_resolution=0.8` | `cell_types` (the same function types RNA and ATAC cells) |
| `rna_sample_level_batch_col: "batch"` | `batch_col=["batch"]` | `compute_sample_embedding`, and RAISIN / trajectory DGE |
| `multiomics_consistency_threshold: 0.05` | `consistency_threshold=0.05` | `multiomics_preparation` (scGLUE) |

The names are not always parallel across modalities; see [section 7](#7-multi-omics-keys-that-differ). The [Config key reference](config_reference.md) lists the consuming function and parameter for all 304 keys.

## 4. The stage graph

Each pipeline (RNA, ATAC, multi-omics) runs in two phases: per-modality preparation, then the shared downstream analysis. Each stage has an on/off flag. The tree below is the RNA pipeline; ATAC is identical with `atac_` keys, and multi-omics differs only in phase 1. Folders are relative to `{output_dir}/rna/` (or `rna_output_dir` if set).

```text
run_rna_pipeline
│  Phase 1: rna_wrapper
├── rna_preprocessing ................ preprocess                → preprocess/adata_preprocessed.h5ad
├── rna_cell_type_cluster ............ cell_types                → obs['cell_type'], preprocess/cell_type.csv, umap_cell_type.png
├── rna_derive_sample_embedding ...... compute_sample_embedding  → sample_embedding/sample_embedding.csv, sample_embedding_blocks.npz
│     └── rna_autotune_enable ........ run_autotune instead      → sample_embedding/ (+ autotune_record.txt)
│        (no embedding → every downstream step below is skipped)
│  Phase 2: downstream_analysis, in this order
├── rna_sample_distance_calculation .. sample_distance           → Sample_distance/<method>/
├── rna_trajectory_analysis
│     ├── rna_trajectory_supervised: true ... CCA_Call          → CCA/  (+ CCA_test/ if rna_cca_pvalue)
│     ├── rna_trajectory_supervised: false .. TSCAN             → TSCAN/
│     └── rna_trajectory_dge ......... GAM on the pseudotime     → trajectoryDEG/
├── rna_sample_cluster ............... k-means on the embedding  → sample_cluster/
├── rna_proportion_test .............. proportion_test           → sample_cluster/proportion_test/
├── rna_cluster_dge .................. raisinfit + run_pairwise_tests → raisin_results/
├── rna_visualize_data ............... visualization             → visualization/
└── rna_dimension_association_analysis  run_dimension_association_analysis → sample_association/
```

Plus `{output_dir}/sys_log/main_process_status.json`, which records which stages ran. What each stage computes is in the [Quick start](quickstart.md); every file, and how to load it back, is in [Understanding the outputs](outputs.md).

Four dependencies are easy to miss:

1. **Trajectory DGE runs only inside the trajectory step.** `rna_trajectory_dge: true` does nothing when `rna_trajectory_analysis: false`. It fits the pseudotime from CCA or TSCAN, whichever ran.
2. **The proportion test and RAISIN compare either k-means clusters or a metadata column, never both.** If `*_cluster_differential_gene_group_col` is set, its groups **replace** the k-means clusters. The demo sets it to `"sev.level"`, so the demo compares severity groups, not clusters. Set it to `null` to compare the clusters from `*_sample_cluster` (which must then be `true`; with neither, both steps are skipped).
3. **Phase 2 needs a sample embedding.** With `rna_derive_sample_embedding: false`, the embedding must already be stored in the cell-level file (see [section 9](#9-resuming-a-run)); otherwise every downstream stage is skipped with a message.
4. **Some failures stop the run, others are only printed.** An exception in preprocessing, cell typing, the embedding, sample distance, trajectory (including DGE), clustering or visualization stops the run (with `continue_on_error: true`, only that modality). Errors in the proportion test, RAISIN, the dimension association and the multi-omics embedding plots are caught, printed, and the run continues.

The multi-omics phase 1 is scGLUE integration (`multiomics_integration` and the five `multiomics_run_glue_*` sub-flags), a Harmony pass that builds the sample-removed `Z_comp` from `X_glue`, joint cell typing (`multiomics_cell_type_cluster`) and the sample embedding (`multiomics_derive_sample_embedding`). Its phase 2 is the same as above except that `visualize_data` is always off, `multiomics_visualize_embedding` draws the embedding plots instead, and trajectory DGE and RAISIN use the RNA cells only. See [Multi-omics from a config](multiomics_config.md).

## 5. Global keys

| Key | Default | What it really does |
|---|---|---|
| `output_dir` | (required) | Root folder. Each pipeline writes to `{output_dir}/rna`, `/atac` or `/multiomics` unless `rna_output_dir`, `atac_output_dir` or `multiomics_output_dir` is set. |
| `run_rna_pipeline` | **`true`** | Run the RNA pipeline. Note the default: in Python, set it to `false` for an ATAC-only run. |
| `run_atac_pipeline`, `run_multiomics_pipeline` | `false` | Run the ATAC / multi-omics pipelines. Order is always RNA, ATAC, multi-omics. |
| `use_gpu` | `false` | Run the cell-level steps (Harmony, neighbor graph, Leiden, UMAP) on the GPU through RAPIDS. Nothing switches to the GPU unless this is `true`; multi-omics also needs `multiomics_use_gpu: true`. If RAPIDS cannot initialise, SampleDisco warns and runs on the CPU, and `gpu_available` in `sys_log/main_process_status.json` records which backend ran. See [GPU acceleration](gpu.md). |
| `initialization` | **`true`** | **Has no effect on which stages run.** `true` deletes `{output_dir}/{rna,atac,multiomics}/result/`, a folder no current stage writes, and rewrites the status JSON. `false` reloads the status JSON and prints it. Which stages run is decided only by the stage flags. |
| `verbose` | `true` | Progress messages (multi-omics uses `multiomics_verbose`). |
| `save_intermediate` | `true` | Read **only by multi-omics**: re-save the union `adata_sample.h5ad` after the Harmony pass and after the embedding. RNA and ATAC always write their files. |
| `large_data_need_extra_memory` | `false` | **GPU only**: use managed (oversubscribable) GPU memory instead of a memory pool. No effect on CPU. |
| `continue_on_error` | `false` | `false`: the first failing pipeline raises and the run stops. `true`: the error is stored in `results["rna_error"]` (or `atac_error`, `multiomics_error`), the next pipeline still runs, and the CLI exits 1 at the end. Use it for batch jobs where partial results are useful. |
| `random_state` | `42` | Seeds the RNA and ATAC sample embedding and RAISIN, **nothing else**. Preprocessing, Harmony and cell typing use a fixed seed (42), sample k-means clustering seed 0, and multi-omics its own `multiomics_random_state` (which seeds only the Harmony pass). |

Runs are not bit-for-bit reproducible, even on the same node with the same seeds: six runs of the demo template (three nodes, 8 cores each) gave RNA CCA scores between 0.981 and 0.989 and ATAC scores between 0.913 and 0.993, with different best PC pairs. Compare results across runs at the level of conclusions, not decimals.

## 6. Inputs, columns and batches

=== "RNA / ATAC"

    | Key | Purpose |
    |---|---|
    | `rna_count_data_path` | Cell-level `.h5ad` with **raw counts** in `.X`. Required whenever the pipeline runs, even when resuming. |
    | `rna_sample_meta_path` | Optional sample-level CSV: one row per sample, a column named like `rna_sample_col`, any other columns. They are joined onto every cell; same-named `.obs` columns are replaced. |
    | `rna_cell_meta_path` | Optional cell-level CSV with a `barcode` column matching `obs_names`. |
    | `rna_output_dir` | Overrides `{output_dir}/rna`. |
    | `rna_adata_cell_path` | Resume input, see [section 9](#9-resuming-a-run). |

    The ATAC keys have the same names with `atac_`. Leave the metadata paths `null` when the columns are already in `.obs`, as in the demo data. The file formats are described in [Using your own data](own_data.md).

=== "Multi-omics"

    | Key | Purpose |
    |---|---|
    | `multiomics_rna_file`, `multiomics_atac_file` | Raw RNA and ATAC `.h5ad` files. **Both must be set** whenever the multi-omics pipeline runs, even with `multiomics_integration: false`. |
    | `multiomics_rna_sample_meta_file`, `multiomics_atac_sample_meta_file` | Optional sample-level CSVs, one per modality. Note the `_file` suffix, not `_path`. |
    | `multiomics_additional_hvg_file` | Optional text file of extra genes to keep through the scGLUE preprocessing. |
    | `multiomics_integrated_h5ad_path` | An already integrated object (cells of both modalities, with `X_glue`), read when `{multiomics output}/preprocess/adata_sample.h5ad` does not exist yet. It is input only: every update goes to that output file. |
    | `multiomics_output_dir` | Overrides `{output_dir}/multiomics`. |

**Columns.**

- `*_sample_col` is the per-cell sample ID. A column literally named `sample` is the safest choice (see [Using your own data](own_data.md)).
- `*_celltype_col` names the column that holds cell-type labels. **Cell typing always writes its labels to `obs['cell_type']`**, whatever this key says, and once it has run, the embedding and every downstream step use `obs['cell_type']`. The key matters in two cases: with `*_existing_cell_types: true` it is the column your labels are copied from, and with `*_cell_type_cluster: false` it is the column the embedding reads. In 0.5.0 the embedding always read `*_celltype_col`, so labels merged by `*_n_target_cell_clusters` were ignored unless the column was `cell_type`.
- `multiomics_modality_col` separates RNA and ATAC cells in the integrated object. Joint cell typing expects the values `RNA` and `ATAC`.

**Batches.** There are two kinds of batch key, and only RNA and ATAC have both.

| Key | Type | Where it is used |
|---|---|---|
| `rna_cell_level_batch_key`, `atac_cell_level_batch_key` | list of columns | Cell-level Harmony in preprocessing. The sample column is added automatically for `Z_comp`. With no non-sample column here, `Z_rmd` is the uncorrected PCA (LSI for ATAC). |
| `rna_sample_level_batch_col`, `atac_sample_level_batch_col` | a column **or a list of columns** | Sample-level Harmony on the final embedding, RAISIN's ComBat (first column only) and the trajectory-DGE pseudobulk. With two or more columns, for example `rna_sample_level_batch_col: ["batch", "site"]`, Harmony corrects all of them jointly (multi-covariate Harmony). |
| `multiomics_batch_col` | one column | scGLUE's batch, the `Z_comp` Harmony pass and sample-level Harmony. There is no multi-omics cell-level batch key. |

Sample-level Harmony runs only when there are **at least 8 units** (samples, or sample × modality units in multi-omics) **and more than one batch value**; otherwise the embedding is left uncorrected without a warning. `*_sample_embedding_batch_method` chooses `"harmony"` (default), `"linear"` (regress the batch out) or `"none"`. [Sample embedding in depth](sample_embedding.md) explains the blocks, weights and batch correction.

## 7. Multi-omics keys that differ

Most multi-omics keys mirror the RNA ones with a `multiomics_` prefix. These do not:

| RNA / ATAC key | Multi-omics equivalent |
|---|---|
| `*_sample_meta_path` | `multiomics_rna_sample_meta_file`, `multiomics_atac_sample_meta_file` |
| `*_sample_level_batch_col` (string or list) | `multiomics_batch_col` (one column) |
| `*_cell_level_batch_key` | none |
| `*_leiden_cluster_resolution` | `multiomics_cluster_resolution` |
| `*_umap` | `multiomics_generate_umap_celltype` |
| `*_cell_embedding_column` | `multiomics_use_rep_celltype` (cell typing only) |
| `*_n_target_cell_clusters` | none: `multiomics_n_target_clusters` exists but is **ignored** |
| `*_existing_cell_types` | none: `multiomics_existing_cell_types` is **ignored**. To keep labels, set `multiomics_cell_type_cluster: false` with labels already in `multiomics_celltype_col`. |
| `*_min_cells`, `*_min_genes`, … | `multiomics_rna_min_cells`, `multiomics_atac_min_features`, … (QC of the per-modality files used by DGE and RAISIN) |
| `*_visualize_data`, `*_plot_*_flag`, `*_age_bin_size`, `*_age_column` | none: the `multiomics_` versions are **ignored** (multi-omics never runs `visualization`). Use `multiomics_visualize_embedding` and `multiomics_visualization_grouping_column`; without the latter the plots are uncoloured. |
| `verbose`, `random_state`, `use_gpu` | `multiomics_verbose`, `multiomics_random_state`, `multiomics_use_gpu` (GPU needs **both** `use_gpu` and `multiomics_use_gpu`) |

The trajectory keys (`*_trajectory_col`, `*_trajectory_supervised`, `*_n_cca_pcs`, `*_cca_pvalue`, `*_tscan_n_clusters`, `*_tscan_origin`, `*_tscan_pseudotime_mode`, `*_trajectory_visualization_label`) exist for all three pipelines. `*_tscan_origin` is an integer, the index of the endpoint cluster TSCAN starts from (`null` picks one). `*_tscan_pseudotime_mode` is `"rank"` or `"distance"`; anything else raises.

## 8. A worked change: unsupervised trajectory, cluster-based tests

The demo config fits a supervised trajectory (CCA against `sev.level`) and runs the proportion test on the severity groups. This change switches both to the unsupervised versions, adds two composition-based sample distances and turns ATAC off. Edit these keys in `my_config.yaml` (a fresh copy from step 1) and leave the rest as it is:

```yaml
output_dir: "./covid_rna_tscan"
run_atac_pipeline: false

rna_trajectory_supervised: false        # TSCAN instead of CCA
rna_tscan_n_clusters: 3
rna_tscan_pseudotime_mode: "distance"   # or "rank"

rna_cluster_differential_gene_group_col: null   # proportion test on the k-means clusters
rna_cluster_number: 3

rna_sample_distance_methods: ["cosine", "chi_square", "jensen_shannon"]
```

`*_sample_distance_methods` accepts the SciPy metrics `euclidean`, `sqeuclidean`, `minkowski`, `cityblock`, `chebyshev`, `cosine`, `correlation`, `hamming`, `jaccard`, `canberra`, `braycurtis` and `matching`, computed on the sample embedding, and `EMD`, `chi_square` and `jensen_shannon`, computed from each sample's cell-type composition. An unknown name is skipped with a warning. In 0.5.0, `chi_square` and `jensen_shannon` failed with a `NameError`.

```bash
sampledisco -m complex --config my_config.yaml
```

This took 1 min 33 s with a peak of 3.3 GB RAM on 8 CPU cores (2 min 32 s on a slower node).

Writes → `covid_rna_tscan/rna/` with `TSCAN/` in place of `CCA/`: `TSCAN/X_DR_sample_pseudotime.csv`, `clusters_by_cluster_X_DR_sample.png` and `clusters_by_grouping_X_DR_sample.png` (coloured by `rna_trajectory_visualization_label`). The two composition distances go to `Sample_distance/chi_square/proportion_DR_distance/` and `Sample_distance/jensen_shannon/proportion_DR_distance/` (distance matrix, heatmap, `cell_type_abundances.pdf` and the group-separation score). The other folders are as in [section 4](#4-the-stage-graph).

- TSCAN found 3 clusters of 5, 2 and 1 samples and ordered the healthy donors first: four of the five lowest pseudotimes are `HD-*` samples, the three highest `CoV-*` samples.
- Trajectory DGE fitted 9,829 of 15,685 cell type × gene features on that pseudotime, none significant (the GAM is underpowered on 8 samples and says so in a warning).
- The group-separation score (mean between-group / within-group distance for `sev.level`) was 1.42 for cosine, 1.15 for chi-square and 1.09 for Jensen-Shannon.
- The k-means clusters had 6, 1 and 1 samples, so the proportion test skipped all three comparisons (`proportion_test_significant_summary.txt` lists them under `Skipped comparisons (degenerate group sizes …)`) and wrote only the heatmaps. With 8 samples, cluster-based tests rarely have two samples per group; on a larger cohort they do.

## 9. Resuming a run

There is no "resume from where it stopped" switch: `initialization` and the status file do not skip anything. You resume by turning finished stages **off**, and the pipeline then reads their results from disk. The rules, per modality:

| Stage turned off | What the pipeline loads instead |
|---|---|
| `rna_preprocessing: false` | The cell-level file `rna_adata_cell_path`, or by default `{rna output}/preprocess/adata_preprocessed.h5ad`. A missing file raises `ValueError: Preprocessed data not found at …`. |
| `rna_cell_type_cluster: false` | The labels already in that file, read from `rna_celltype_col`. Leave it `"cell_type"` to reuse the labels an earlier run wrote. |
| `rna_derive_sample_embedding: false` | The embedding stored in that file's `uns['X_DR_sample']`. It is there only if an earlier run had `rna_save_cell_adata_after_embedding: true` (the default). Without it, every downstream stage is skipped. |
| downstream flags | Each downstream stage reruns from the embedding when its flag is `true`; nothing downstream is loaded from disk. |
| `multiomics_integration: false` | The union `{multiomics output}/preprocess/adata_sample.h5ad` if an earlier run wrote it, otherwise `multiomics_integrated_h5ad_path`. Every update (Harmony pass, cell types, embedding) is written to the output file, never to the input. The embedding is read from `uns['X_DR_sample']` when `multiomics_derive_sample_embedding: false`. Because the output file wins, pointing `multiomics_integrated_h5ad_path` at a different object has no effect until you change `output_dir` or delete that file. |

!!! warning "Resume keys that do nothing"
    `rna_adata_sample_path`, `rna_pseudo_adata_path`, `atac_adata_sample_path`, `atac_pseudo_adata_path` and `multiomics_pseudobulk_h5ad_path` are accepted, required by the CLI, and **never read**. Earlier versions of this page told you to set them; leave them `null`. Ten more keys are equally dead (all multi-omics): they are flagged in the [Config key reference](config_reference.md#config-key-reference).

For example, to add RAISIN and a supervised trajectory to the run of the previous section without recomputing preprocessing, cell types or the embedding, change these keys in the same `my_config.yaml`:

```yaml
rna_preprocessing: false              # load covid_rna_tscan/rna/preprocess/adata_preprocessed.h5ad
rna_cell_type_cluster: false          # keep the cell types stored in that file
rna_derive_sample_embedding: false    # reuse its uns['X_DR_sample']
rna_adata_cell_path: null             # null = the default location above

rna_trajectory_supervised: true       # now CCA against sev.level
rna_cluster_dge: true                 # and add RAISIN
```

```bash
sampledisco -m complex --config my_config.yaml
```

The log shows `Resuming from previous progress:` and goes straight from `Starting RNA wrapper function...` to the downstream steps; nothing under `preprocess/` or `sample_embedding/` is rewritten. The distance scores are identical to the previous run (same embedding), CCA scored 0.998 on PC1 + PC3, and RAISIN compared the three k-means clusters (it warns that the variance of the two single-sample clusters cannot be estimated).

Writes → new `CCA/` and `raisin_results/` folders (`raisin_results/<a>_vs_<b>/raisin_results.csv` and volcano plots for each pair, plus `summary_plots/`), and fresh copies of the distance, trajectory-DGE, clustering, visualization and association outputs. The older `TSCAN/` folder stays. It took 1 min 20 s (2 min 22 s on a slower node), but RAISIN raised the peak memory from 3.3 GB to **15.0 GB**; see [Using your own data](own_data.md#8-memory-and-runtime).

!!! note "Things that change files in place"
    - A resumed run that keeps `*_cell_type_cluster: true` re-types the cells and **overwrites** `preprocess/adata_preprocessed.h5ad`. If `*_derive_sample_embedding` is `false` at the same time, the stored embedding was computed from the old labels.
    - The embedding step updates `uns['X_DR_sample']` inside the cell-level file (`*_save_cell_adata_after_embedding: false` keeps only the CSV).
    - Multi-omics writes its cell types and embedding to `{multiomics output}/preprocess/adata_sample.h5ad`. In 0.5.0 it wrote them **into the file named by `multiomics_integrated_h5ad_path`** (with the demo config, `data/test_multiomics_integrated.h5ad`); keep a copy of that file if you run 0.5.0.

## 10. Python with a partial config

`wrapper()` has a default for every key except `output_dir`, so from Python you can pass only what you need. Everything you leave out takes the `wrapper()` default, which is the "Default" column of the [Config key reference](config_reference.md), not the demo value.

```python
from sampledisco.wrapper.wrapper import wrapper

res = wrapper(
    output_dir="partial_run",
    rna_count_data_path="data/test_RNA.h5ad",
    rna_cell_level_batch_key=["batch"],      # cell-level Harmony on batch
    rna_sample_level_batch_col="batch",      # sample-level Harmony on batch
    rna_trajectory_supervised=True,          # CCA against rna_trajectory_col ("sev.level")
    rna_n_cca_pcs=2,
    rna_spline_order=2,                      # order 3 fits no gene on 8 samples
)
pseudo = res["rna_results"]["pseudo_adata"]
print(pseudo)
print(pseudo.obs[["sev.level", "batch", "pseudotime_sample"]])
```

```text
AnnData object with n_obs × n_vars = 8 × 7
    obs: 'celltype', 'sample', 'old', 'large', 'batch', 'sev.level', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'pseudotime_sample', 'cluster_sample_kmeans'
    uns: 'X_DR_sample', 'sample_embedding_params'
    obsm: 'X_DR_sample'
              sev.level  batch  pseudotime_sample
sample                                           
CoV-3-Wilk            4   Wilk           0.466959
CoV-34-Aruna          4  Aruna           0.000000
CoV-36-Aruna          4  Aruna           0.660007
CoV-4-Wilk            4   Wilk           0.379138
HD-1-Wilk             1   Wilk           0.740808
HD-20-Aruna           1  Aruna           0.702005
HD-21-Aruna           1  Aruna           0.758945
HD-4-Wilk             1   Wilk           1.000000
```

Writes → `partial_run/rna/` with `preprocess/`, `sample_embedding/`, `Sample_distance/`, `CCA/`, `trajectoryDEG/`, `sample_cluster/` and `visualization/`, plus `partial_run/sys_log/`. It took 1 min 26 s with a peak of 3.6 GB RAM on 8 CPU cores (2 min 36 s on a slower node).

The defaults that differ from the demo template show up in the result:

- RNA only, because `run_atac_pipeline` defaults to `false` (but `run_rna_pipeline` to `true`).
- `rna_min_cells=500` and `rna_min_genes=500`: preprocessing kept 29,680 cells and 9,936 genes (the template's 100 / 200 keep 29,988 cells and 12,301 genes).
- Leiden at resolution 0.8 without a target count gave 19 cell types (the template asks for 8).
- `rna_cluster_number=4` k-means clusters (`cluster_sample_kmeans`); no proportion test, RAISIN or association analysis, because those flags default to `false`.
- The grouping and trajectory columns silently fell back to `sev.level`.
- With CCA on 2 PCs (score 0.75; the template's 5 PCs give 0.98–0.99), the pseudotime does not order the two severity groups cleanly. The 8-sample demo is too small to judge a trajectory; see the [Sound Life RNA tutorial](soundlife_rna.md) for a 32-donor cohort.

## 11. Pitfalls

!!! warning "The summary CSV is deleted"
    `*_summary_sample_csv_path` names a CSV that the sample-distance step appends its group-separation scores to. The pipeline deletes that file again right after the trajectory step, so after a run it is gone.

    - **0.6.0** always deletes the default `{modality output}/summary_sample.csv`. A path you supply is deleted only if the run created it; a file that existed before the run is kept, with the scores appended.
    - **In 0.5.0** the file is deleted in every case, including one that existed before the run.

    The per-method scores stay in `Sample_distance/<method>/distance_check_results_*.txt` either way.

!!! warning "Small cohorts"
    The 8-sample demo is at the lower limit of several steps: the embedding has at most n_samples − 1 PCs, sample-level Harmony needs at least 8 units, the CCA score is close to 1 by construction when `*_n_cca_pcs` approaches the number of samples (keep `*_cca_pvalue: false` on the demo), and trajectory DGE needs `*_spline_order: 2`. [Using your own data](own_data.md#6-small-cohorts) has the exact limits.

!!! warning "GPU and CPU runs differ"
    `use_gpu: true` does not give the same numbers as a CPU run. The sample embedding's k-means runs on scikit-learn on both backends (cuML 24.12 has no `MiniBatchKMeans`); the differences come from GPU cell typing (RAPIDS neighbors, Leiden and PCA, which are not deterministic from run to run) and from Harmony running on a CUDA device. Pick one backend per study. [GPU acceleration](gpu.md) has the measured differences and run times.

!!! note "Autotune ignores some keys"
    With `*_autotune_enable: true` the embedding is computed by `run_autotune` (CPU only), which ignores `*_sample_embedding_use_clr`, `*_use_rmd`, `*_block_weights` and `*_rmd_weight`. Without a sample-level batch column and without `*_autotune_grouping_col`, it does not search at all. See [Autotune](downstream/autotune.md).

!!! note "Example configs in the repository"
    The SampleDisco repository has more configs under `config/` (for example `config/config_covid_rna.yaml`). They validate against 0.5.0 and 0.6.0, but their paths point to the authors' cluster.
