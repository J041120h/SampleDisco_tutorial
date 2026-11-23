# Configuration guide

A SampleDisco run is one YAML file. The CLI reads it into a flat dictionary and calls [`wrapper(**config)`](../api/pipeline/wrapper.md), so every key is a parameter of `wrapper()`. This page explains how the keys fit together. Every key, with its default and the function that reads it, is in the [Config key reference](config_reference.md).

!!! warning "Change these keys before you run on your own data"
    The template is wired to the COVID-19 demo, whose phenotype column is `sev.level`. These keys (each has an `rna_`, `atac_` and `multiomics_` version) name that column, and **`null` does not switch them off: it falls back to `sev.level`.**

    | Key | Demo value | If the column is missing |
    |---|---|---|
    | `*_grouping_columns` | `["sev.level"]` | The distance step skips its group score; then the visualization step raises `ValueError: Missing required columns in pseudobulk AnnData` and the run stops. |
    | `*_trajectory_col` | `"sev.level"` | CCA (`*_trajectory_supervised: true`) raises `ValueError: Trajectory column 'sev.level' not found`. For CCA the column must be **numeric**: text values become NaN and those samples are dropped. |
    | `*_trajectory_visualization_label` | `["sev.level"]` | TSCAN skips its group plot, with a warning. |
    | `*_cluster_differential_gene_group_col` | `"sev.level"` | The proportion test prints the `KeyError` and writes nothing; RAISIN falls back to the k-means clusters. |
    | `*_association_categorical_cols`, `*_association_continuous_cols` | `null` (every sample-level column) | Nothing breaks; list the covariates you care about. |

    Also check the data paths, the batch keys and the QC thresholds (`*_min_cells`, `*_min_genes`, …), which are tuned for the demo. [Using your own data](own_data.md) walks through a real example.

## Get a template and run it

```bash
sampledisco --init-config my_config.yaml
```

- This writes the demo config that ships with the package, the same file as [`config_demo.yaml`](../assets/config_demo.yaml). Without a path it writes `./config.yaml`; it never overwrites a file.
- The template is complete and runs as written on the [demo data](demo_data.md) in `./data/`, writing to `./sampledisco_demo_output/`.
- Edit values, but don't delete keys.

=== "CLI"

    ```bash
    sampledisco -m complex --config my_config.yaml
    ```

    - `-m complex` is the only mode and the default.
    - Before anything runs, the CLI checks the file against `wrapper()`. An unknown key raises `Unexpected parameter in config: '<key>'`, a missing one `Missing required parameter in config: '<key>'`. Only the three `*_save_cell_adata_after_embedding` keys may be left out. The values of `*_sample_embedding_rmd_weight` and `*_sample_embedding_batch_method` are checked too.
    - The exit status is 1 if any pipeline failed.

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

    - `results["rna_results"]["pseudo_adata"]` is the sample-level AnnData and `results["rna_results"]["adata"]` the cell-level one. The same goes for `atac_results` and `multiomics_results`.
    - `wrapper()` does not check the key set. A misspelled key raises `TypeError`, but a **missing key silently takes the `wrapper()` default**, which often differs from the demo value (see [below](#python-with-a-partial-config)).

On the demo template (RNA and ATAC) a run takes about 3–15 minutes on 8 CPU cores, depending on the machine, and peaks at about 7.5 GB of RAM.

## How keys map to functions

The prefix says which pipeline reads a key: `rna_`, `atac_` or `multiomics_`, and no prefix for global keys. The rest of the name is usually the parameter of the function the value ends up in: `rna_min_genes: 200` becomes `min_genes=200` in the RNA `preprocess`, and `rna_sample_level_batch_col: "batch"` becomes `batch_col=["batch"]` in `compute_sample_embedding`. The reference lists the function and parameter for every key.

## The stage graph

Each pipeline runs in two phases: per-modality preparation, then the shared downstream analysis. Each stage has an on/off flag. The tree shows the RNA pipeline. ATAC is identical with `atac_` keys. Folders are relative to `{output_dir}/rna/`.

```text
run_rna_pipeline
│  Phase 1: rna_wrapper
├── rna_preprocessing ................ preprocess                → preprocess/adata_preprocessed.h5ad
├── rna_cell_type_cluster ............ cell_types                → obs['cell_type'], preprocess/cell_type.csv
├── rna_derive_sample_embedding ...... compute_sample_embedding  → sample_embedding/
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
├── rna_cluster_dge .................. RAISIN                    → raisin_results/
├── rna_visualize_data ............... visualization             → visualization/
└── rna_dimension_association_analysis  run_dimension_association_analysis → sample_association/
```

Every file, and how to load it back, is in [Outputs and resuming runs](outputs.md). Four dependencies are easy to miss:

1. **Trajectory DGE runs only inside the trajectory step.** `rna_trajectory_dge: true` does nothing when `rna_trajectory_analysis: false`.
2. **The proportion test and RAISIN compare a metadata column or the k-means clusters, never both.** If `*_cluster_differential_gene_group_col` is set (the demo sets `"sev.level"`), its groups replace the clusters. Set it to `null` to compare the clusters from `*_sample_cluster`, which must then be `true`.
3. **Phase 2 needs a sample embedding.** With `*_derive_sample_embedding: false` it is read from the cell-level file ([resuming](outputs.md#resuming-a-run)); without one, every downstream stage is skipped with a message.
4. **Some failures stop the run, others are only printed.** An error in preprocessing, cell typing, the embedding, distances, trajectory, clustering or visualization stops the run (with `continue_on_error: true`, only that modality). Errors in the proportion test, RAISIN, the association analysis and the multi-omics plots are printed, and the run continues.

## Global keys

| Key | Default | What it does |
|---|---|---|
| `output_dir` | (required) | Root folder. Each pipeline writes to `{output_dir}/rna`, `/atac` or `/multiomics`, unless `rna_output_dir`, `atac_output_dir` or `multiomics_output_dir` is set. |
| `run_rna_pipeline`, `run_atac_pipeline`, `run_multiomics_pipeline` | `true`, `false`, `false` | Which pipelines run, always in this order. From Python, pass `run_rna_pipeline=False` for an ATAC-only run. |
| `use_gpu` | `false` | Runs the cell-level steps on the GPU; multi-omics also needs `multiomics_use_gpu`. See [GPU acceleration](gpu.md). |
| `large_data_need_extra_memory` | `false` | GPU only: CUDA managed memory instead of a memory pool. |
| `initialization` | `true` | Resets (`true`) or reloads (`false`) the [status file](outputs.md#the-status-file). It never makes a stage run or skip. |
| `continue_on_error` | `false` | `true`: a failing pipeline stores its error in `results["rna_error"]` (or `atac_error`, `multiomics_error`) and the next pipeline still runs. The CLI still exits 1. |
| `random_state` | `42` | Seeds the RNA and ATAC sample embedding and RAISIN only. Preprocessing and cell typing use a fixed seed; multi-omics uses `multiomics_random_state`. |
| `save_intermediate` | `true` | Read only by multi-omics: re-saves its cell-level file after the Harmony pass and after the embedding. |

Results can differ slightly between runs and machines; see [run-to-run variation](own_data.md#run-to-run-variation).

## Inputs, columns and batches

**Inputs.**

- `rna_count_data_path`, `atac_count_data_path`: the cell-level `.h5ad` with **raw counts** in `.X`. Required whenever that pipeline runs, even when resuming.
- `*_sample_meta_path`, `*_cell_meta_path`: optional CSVs of sample-level and cell-level metadata. The formats are in [Using your own data](own_data.md#sample-and-cell-metadata-files). Leave them `null` when the columns are already in `.obs`.
- `*_adata_cell_path`: the input for a resumed run ([resuming](outputs.md#resuming-a-run)).

**Columns.**

- `*_sample_col` is the per-cell sample ID. A column literally named `sample` is the safest choice.
- `*_celltype_col`: cell typing always writes its labels to `obs['cell_type']`, and the embedding and every downstream step then use that column. This key matters in two cases: with `*_existing_cell_types: true` it is the column your labels are copied from, and with `*_cell_type_cluster: false` it is the column the embedding reads.

**Batches.** There are two kinds of batch key:

| Key | Type | Where it is used |
|---|---|---|
| `rna_cell_level_batch_key`, `atac_cell_level_batch_key` | list of columns | Cell-level Harmony in preprocessing. The sample column is added automatically for `Z_comp`. |
| `rna_sample_level_batch_col`, `atac_sample_level_batch_col` | a column or a list | Sample-level Harmony on the embedding (two or more columns: multi-covariate Harmony), RAISIN's ComBat (first column) and the trajectory-DGE pseudobulk. |
| `multiomics_batch_col` | one column | scGLUE's batch, the `Z_comp` Harmony pass and the sample-level Harmony. |

Sample-level Harmony runs only with **at least 8 units and more than one batch value**; otherwise the embedding is left uncorrected. `*_sample_embedding_batch_method` chooses `"harmony"` (default), `"linear"` or `"none"`. [Sample embedding in depth](sample_embedding.md#pca-and-sample-level-batch-correction) explains when each one helps.

## Multi-omics

The multi-omics pipeline (`run_multiomics_pipeline: true`) has its own phase 1: scGLUE integration, a Harmony pass that turns `X_glue` into the sample-removed `Z_comp`, joint cell typing, and a sample embedding with one unit per sample × modality. Phase 2 is the stage graph above, with two changes: `multiomics_visualize_embedding` draws the embedding plots in place of the visualization step, and trajectory DGE and RAISIN use only the RNA cells. Joint cell typing expects the values `RNA` and `ATAC` in `multiomics_modality_col`.

**Start from a pre-integrated file.** A fresh template (`sampledisco --init-config`) already has `multiomics_integration: false` and points `multiomics_integrated_h5ad_path` at the demo's scGLUE output. Switch the pipelines and colour the plots:

```yaml
run_rna_pipeline: false
run_atac_pipeline: false
run_multiomics_pipeline: true
multiomics_integration: false
multiomics_integrated_h5ad_path: "./data/test_multiomics_integrated.h5ad"
multiomics_visualization_grouping_column: ["sev.level"]
```

- `multiomics_integration: false` skips scGLUE and reads `multiomics_integrated_h5ad_path`: the cells of both modalities with their `X_glue` embedding. scGLUE, `multiomics_rna_file` and `multiomics_atac_file` are then not needed. Without `multiomics_visualization_grouping_column` the plots are uncoloured. On the demo the run takes 1–5.5 minutes on 8 CPU cores and under 3 GB of RAM.
- **The input is never modified.** Every update goes to `{multiomics output}/preprocess/adata_sample.h5ad`. The pipeline reads that file first when it exists, so a second run continues from the first. Delete it, or change `output_dir`, to start again from the input.
- **DGE and RAISIN are skipped.** They need the per-modality RNA expression file `preprocess/adata_rna_preprocessed.h5ad`, which only an integration run writes; the integrated object holds embeddings, not expression. With `multiomics_trajectory_dge` or `multiomics_cluster_dge` on, the log says `WARNING: skipping multi-omics trajectory DGE and cluster DGE (RAISIN)` and every other step runs.

**Train scGLUE from the config.** Set `multiomics_integration: true`, leave `multiomics_integrated_h5ad_path` empty, point `multiomics_rna_file` and `multiomics_atac_file` at the raw counts, and preferably give the run its own `multiomics_output_dir`. This needs scGLUE and `bedtools` (see the [multi-omics tutorial](multiomics.md)) and writes the per-modality files, so DGE and RAISIN can run. On the demo it took about 30 minutes on a GPU node and 15 GB of RAM. The scGLUE keys:

| Key | Template | Meaning |
|---|---|---|
| `multiomics_rna_file`, `multiomics_atac_file` | `./data/test_RNA.h5ad`, `./data/test_ATAC.h5ad` | raw counts; `.obs` must have the sample column |
| `multiomics_ensembl_release`, `multiomics_species` | `98`, `homo_sapiens` | gene annotation for the guidance graph, downloaded on first use |
| `multiomics_n_top_genes`, `multiomics_n_top_peaks` | `2000`, `50000` | highly variable genes and peaks given to scGLUE |
| `multiomics_glue_max_epochs`, `multiomics_glue_data_batch_size` | `50`, `128` | training length and minibatch size |
| `multiomics_run_glue_twice_for_sample_removal` | `false` | `true` trains scGLUE a second time to get `Z_comp`, instead of the Harmony pass |
| `multiomics_batch_col` | empty | a technical batch; scGLUE removes it only when both files have the column |

**GPU.** `multiomics_use_gpu: true` alone does nothing. Multi-omics uses the GPU only when `use_gpu` is also `true`. The template has `use_gpu: false`, so the log says `GPU mode: Disabled`. scGLUE training is the exception: it uses a visible CUDA GPU on its own, whatever these keys say.

**Keys that differ from RNA and ATAC.** Most multi-omics keys mirror the RNA ones with a `multiomics_` prefix. These don't:

| RNA / ATAC key | Multi-omics |
|---|---|
| `*_sample_meta_path` | `multiomics_rna_sample_meta_file`, `multiomics_atac_sample_meta_file` |
| `*_sample_level_batch_col` (column or list) | `multiomics_batch_col` (one column) |
| `*_cell_level_batch_key` | none |
| `*_leiden_cluster_resolution`, `*_umap` | `multiomics_cluster_resolution`, `multiomics_generate_umap_celltype` |
| `*_existing_cell_types`, `*_n_target_cell_clusters` | none: the `multiomics_` versions are ignored. To keep your labels, set `multiomics_cell_type_cluster: false` with the labels in `multiomics_celltype_col`. |
| `*_visualize_data` and the plot flags | none: use `multiomics_visualize_embedding` and `multiomics_visualization_grouping_column` (without it the plots are uncoloured). `multiomics_target_modality` (template `ATAC`) adds plots that highlight one modality. |
| `*_min_cells`, `*_min_genes`, … | `multiomics_rna_min_cells`, `multiomics_atac_min_features`, …: QC of the per-modality files |
| `verbose`, `random_state` | `multiomics_verbose`, `multiomics_random_state` |

From Python, the `cell_types_multiomics` arguments `cluster_resolution`, `transfer_metric` and `use_rep` map to `multiomics_cluster_resolution`, `multiomics_metric` and `multiomics_use_rep_celltype`; `k_neighbors` is fixed at 3. The `compute_sample_embedding` options are the `multiomics_sample_embedding_*` keys.

## A worked change: TSCAN and cluster-based tests

The demo fits a supervised trajectory (CCA against `sev.level`) and runs the proportion test on the severity groups. To switch both to their unsupervised versions and add two composition-based distances, edit these keys in a fresh `my_config.yaml`:

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

```bash
sampledisco -m complex --config my_config.yaml
```

This takes about 2 minutes. `covid_rna_tscan/rna/` then has `TSCAN/` in place of `CCA/` (`X_DR_sample_pseudotime.csv` and two cluster plots), and the composition distances in `Sample_distance/chi_square/` and `Sample_distance/jensen_shannon/`. With 8 samples, k-means gives clusters of one sample, and the proportion test skips every comparison involving them.

`*_sample_distance_methods` accepts SciPy metrics computed on the embedding (`euclidean`, `cosine`, `correlation`, `cityblock`, … 12 in all) and `EMD`, `chi_square` and `jensen_shannon`, computed from each sample's cell-type composition. Unknown names are skipped with a warning.

## Python with a partial config

`wrapper()` has a default for every key except `output_dir`, so from Python you can pass only the keys you need. Everything else takes the `wrapper()` default (the "Default" column of the [Config key reference](config_reference.md)), not the demo value. The differences that matter most: `run_rna_pipeline` defaults to `true` and the other two pipelines to `false`; `rna_min_cells` and `rna_min_genes` default to 500 (the template uses 100 and 200); `rna_n_target_cell_clusters` is `null`, so Leiden decides the number of cell types; `rna_trajectory_supervised` is `false` (TSCAN); the proportion test, RAISIN and the association analysis are off; and the `sev.level` keys in the box at the top fall back to `sev.level` silently.

## Other pitfalls

- **Resuming a run** works by switching finished stages off. The rules and the keys that look like resume switches but are never read are in [Outputs and resuming runs](outputs.md#resuming-a-run).
- **`*_summary_sample_csv_path`** names a working file of the distance step that is deleted at the end of the run. The scores stay in `Sample_distance/<method>/`.
- **Autotune** (`*_autotune_enable: true`) replaces `compute_sample_embedding` and ignores `*_sample_embedding_use_clr`, `*_use_rmd`, `*_block_weights` and `*_rmd_weight`; see [Downstream options and caveats](downstream_details.md).
- **GPU and CPU runs give different numbers**; pick one backend per study ([GPU acceleration](gpu.md)).
- **Small cohorts** hit hard limits in several steps; see [Using your own data](own_data.md#small-cohorts).
