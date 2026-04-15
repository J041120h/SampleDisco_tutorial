# Configuration guide

SampleDisc is driven by a single YAML configuration file. At runtime the file is parsed into a flat dictionary and handed to [`wrapper(**config)`](../api/index.md) — every key in the YAML matches a parameter name of a wrapper function.

The canonical reference is [`code/config/config_covid_rna.yaml`](https://github.com/). This page walks through the same file, block by block, so you can build your own.

## Run it

=== "CLI"

    ```bash
    python /users/hjiang/GenoDistance/code/SampleDisc.py -m complex \
      --config /path/to/your_config.yaml
    ```

=== "Python"

    ```python
    import yaml
    from wrapper.wrapper import wrapper

    with open("/path/to/your_config.yaml") as f:
        config = yaml.safe_load(f)

    wrapper(**config)
    ```

## How the config is organized

A SampleDisc config is a flat YAML file with six logical blocks. Read it top-to-bottom; every key follows the same naming rule:

> **`rna_*`** keys feed the RNA wrapper, **`atac_*`** keys feed the ATAC wrapper, **`multiomics_*`** keys feed the multi-omics wrapper. Anything without a prefix is a global setting.

Example mapping:

| Config key | Becomes | In function |
| --- | --- | --- |
| `rna_preprocessing: true` | `preprocessing=True` | `rna_wrapper(...)` |
| `atac_leiden_cluster_resolution: 0.8` | `leiden_cluster_resolution=0.8` | ATAC `cell_types_linux(...)` |
| `multiomics_consistency_threshold: 0.05` | `consistency_threshold=0.05` | `multiomics_preparation(...)` |

---

## Block A — Global settings

These are top-level keys that apply to the whole run.

| Key | Type | Default | Effect |
| --- | --- | --- | --- |
| `output_dir` | str | — (required) | Root directory. Per-modality outputs land under `{output_dir}/{rna,atac,multiomics}/`. |
| `run_rna_pipeline` | bool | `false` | Run the RNA wrapper branch. |
| `run_atac_pipeline` | bool | `false` | Run the ATAC wrapper branch. |
| `run_multiomics_pipeline` | bool | `false` | Run the multi-omics wrapper branch. |
| `use_gpu` | bool | `false` | Use the `_linux`/GPU-accelerated implementations. Requires RAPIDS + scanpy CUDA wheels. |
| `initialization` | bool | `false` | If `true`, clears existing result folders before running. Use only when you want a clean slate. |
| `verbose` | bool | `true` | Print progress messages. |
| `save_intermediate` | bool | `true` | Write intermediate `.h5ad` outputs to disk. |
| `large_data_need_extra_memory` | bool | `false` | Switch to memory-safe paths for very large datasets. |

```yaml
output_dir: "/dcs07/hongkai/data/harry/result/test"
run_rna_pipeline: true
run_atac_pipeline: false
run_multiomics_pipeline: false
use_gpu: true
initialization: false
verbose: true
save_intermediate: true
```

## Block B — Per-modality inputs and outputs

Each pipeline reads its inputs from dedicated path keys. The expected format is `.h5ad` for counts and `.csv` for metadata.

=== "RNA"

    | Key | Purpose |
    | --- | --- |
    | `rna_count_data_path` | Cell-level `.h5ad` (raw counts preferred). |
    | `rna_sample_meta_path` | Per-sample metadata CSV keyed by `sample`. |
    | `rna_cell_meta_path` | *Optional* per-cell metadata CSV. |
    | `rna_output_dir` | Overrides `{output_dir}/rna`. Leave `null` to use the default. |

=== "ATAC"

    | Key | Purpose |
    | --- | --- |
    | `atac_count_data_path` | Cell-level peak-count `.h5ad`. |
    | `atac_sample_meta_path` | Per-sample metadata CSV. |
    | `atac_cell_meta_path` | *Optional* per-cell metadata CSV. |
    | `atac_output_dir` | Overrides `{output_dir}/atac`. |

=== "Multi-omics"

    | Key | Purpose |
    | --- | --- |
    | `multiomics_rna_file` | RNA `.h5ad`. |
    | `multiomics_atac_file` | ATAC `.h5ad`. |
    | `multiomics_rna_sample_meta_file` / `multiomics_atac_sample_meta_file` | Optional per-modality sample metadata. |
    | `multiomics_additional_hvg_file` | *Optional* plain-text list of extra HVGs to keep through GLUE preprocessing. |
    | `multiomics_output_dir` | Overrides `{output_dir}/multiomics`. |

## Block C — Stage toggles and resume paths

Every pipeline is broken into discrete stages, each guarded by a boolean flag. **Turn a stage off only if you already have its output on disk** — the next stage will read it via a resume path.

### Stage flags

=== "RNA / ATAC"

    ```yaml
    rna_preprocessing: true            # preprocess_linux → adata_cell.h5ad, adata_sample.h5ad
    rna_cell_type_cluster: true        # cell_types_linux  → cell_type column
    rna_derive_sample_embedding: true  # calculate_sample_embedding → pseudobulk_sample.h5ad
    rna_cca_based_cell_resolution_selection: false  # advanced; API-only
    rna_sample_distance_calculation: true
    rna_trajectory_analysis: true
    rna_trajectory_dge: true
    rna_sample_cluster: true
    rna_proportion_test: false
    rna_cluster_dge: false
    rna_visualize_data: true
    ```

=== "Multi-omics"

    ```yaml
    multiomics_integration: true                # multiomics_preparation (GLUE)
    multiomics_integration_preprocessing: true  # integrate_preprocess
    multiomics_cell_type_cluster: true          # cell_types_multiomics_linux
    multiomics_dimensionality_reduction: true   # calculate_multiomics_sample_embedding
    multiomics_find_optimal_resolution: false   # advanced; API-only
    multiomics_sample_distance_calculation: true
    multiomics_trajectory_analysis: true
    multiomics_trajectory_dge: true
    multiomics_sample_cluster: true
    multiomics_proportion_test: false
    multiomics_cluster_dge: false
    multiomics_visualize_embedding: true        # visualize_multimodal_embedding
    # GLUE sub-steps (only relevant while multiomics_integration: true)
    multiomics_run_glue_preprocessing: true
    multiomics_run_glue_training: true
    multiomics_run_glue_gene_activity: true
    multiomics_run_glue_visualization: true
    ```

### Resume paths

If you set a stage flag to `false`, point the next stage at the h5ad it expects:

| Flag turned off | Must provide |
| --- | --- |
| `rna_preprocessing` (or ATAC) | `rna_adata_cell_path`, `rna_adata_sample_path` |
| `rna_derive_sample_embedding` | `rna_pseudo_adata_path` |
| `multiomics_integration` | `multiomics_integrated_h5ad_path` |
| `multiomics_dimensionality_reduction` | `multiomics_pseudobulk_h5ad_path` |

Leave these as `null` when the matching stage is enabled.

## Block D — Common column names

Every modality needs to know which columns in `.obs` identify samples, cell types, and batches.

| Key | Typical value | Meaning |
| --- | --- | --- |
| `*_sample_col` | `"sample"` | Per-cell sample identifier. |
| `*_celltype_col` | `"cell_type"` | Where cell-type assignments are written/read. |
| `*_sample_level_batch_col` | `null` or a string | Sample-level batch for Harmony in sample embedding. |
| `*_cell_level_batch_key` | `null` or a list | Cell-level batch keys for Harmony in preprocessing. |
| `multiomics_modality_col` | `"modality"` | Column that separates RNA and ATAC cells in the integrated object. |

## Block E — Stage-specific tuning knobs

Each stage exposes a small set of numeric/string parameters. The key prefix matches the stage; see the linked API page for full parameter descriptions and defaults.

| Stage | Common keys | API reference |
| --- | --- | --- |
| Preprocess (RNA) | `rna_min_cells`, `rna_min_genes`, `rna_pct_mito_cutoff`, `rna_num_cell_hvgs`, `rna_cell_embedding_num_pcs`, `rna_num_harmony_iterations` | [preprocess_linux](../api/rna/preprocess_linux.md) |
| Preprocess (ATAC) | `atac_min_features`, `atac_max_features`, `atac_num_cell_hvfs`, `atac_cell_embedding_num_pcs`, `atac_tfidf_scale_factor`, `atac_drop_first_lsi`, `atac_doublet_detection` | [preprocess_linux (ATAC)](../api/atac/preprocess_linux.md) |
| Cell types | `*_leiden_cluster_resolution`, `*_n_target_cell_clusters`, `*_existing_cell_types`, `*_umap` | [cell_types_linux](../api/rna/cell_types_linux.md) |
| Sample embedding | `*_sample_hvg_number`, `*_sample_embedding_dimension`, `*_harmony_for_proportion` | [calculate_sample_embedding](../api/shared/calculate_sample_embedding.md) |
| Optimal resolution *(advanced)* | `*_cca_coarse_start/end/step`, `*_cca_fine_range/step`, `*_n_cca_pcs`, `*_cca_compute_corrected_pvalues` | [find_optimal_cell_resolution_linux](../api/shared/find_optimal_cell_resolution_linux.md) |
| Sample distance | `*_sample_distance_methods`, `*_grouping_columns` | [sample_distance](../api/downstream/sample_distance.md) |
| Trajectory | `*_trajectory_supervised`, `*_trajectory_col`, `*_n_cca_pcs`, `*_cca_pvalue`, `*_tscan_origin` | [CCA_Call](../api/downstream/trajectory_cca_call.md) · [TSCAN](../api/downstream/trajectory_tscan.md) |
| Trajectory DGE | `*_fdr_threshold`, `*_effect_size_threshold`, `*_top_n_genes`, `*_num_splines`, `*_spline_order`, `*_trajectory_diff_gene_covariate` | [trajectory GAM DGE](../api/downstream/trajectory_gam_dge.md) |
| Sample clustering | `*_cluster_number` | [cluster](../api/downstream/cluster.md) |
| Proportion test | `*_cluster_differential_gene_group_col` | [proportion_test](../api/downstream/proportion_test.md) |
| Visualization | `*_plot_dendrogram_flag`, `*_plot_cell_type_proportions_pca_flag`, `*_plot_cell_type_expression_umap_flag`, `*_grouping_columns` | [visualization](../api/downstream/visualization.md) |
| Multi-omics integration | `multiomics_ensembl_release`, `multiomics_species`, `multiomics_n_top_genes`, `multiomics_n_pca_comps`, `multiomics_n_lsi_comps`, `multiomics_consistency_threshold`, `multiomics_treat_sample_as_batch`, `multiomics_k_neighbors`, `multiomics_use_rep`, `multiomics_metric` | [multiomics_preparation](../api/multiomics/multiomics_preparation.md) |
| Multi-omics embedding viz | `multiomics_color_col`, `multiomics_visualization_grouping_column`, `multiomics_target_modality`, `multiomics_figsize`, `multiomics_point_size`, `multiomics_alpha`, `multiomics_colormap` | [visualize_multimodal_embedding](../api/downstream/visualize_multimodal_embedding.md) |

## Block F — Output layout

Under `{output_dir}/{modality}/` you will find:

| Directory | Produced by | What's inside |
| --- | --- | --- |
| `preprocess/` | `preprocess_linux` | `adata_cell.h5ad`, `adata_sample.h5ad`, QC summary |
| `pseudobulk/` | `calculate_sample_embedding` | `pseudobulk_sample.h5ad` with `X_DR_expression`, `X_DR_proportion` in `.obsm` |
| `embeddings/` | `calculate_sample_embedding` | CSV exports of the two embeddings |
| `Sample_distance/` | `sample_distance` | Per-metric subdirs with CSVs and heatmap PDFs |
| `CCA/` | `CCA_Call` | 2D CCA plots, contribution plots, pseudotime CSVs |
| `CCA_test/` | `cca_pvalue_test` | Null-distribution plots and p-values |
| `TSCAN/` | `TSCAN` | Cluster-by-group and by-cluster trajectory plots |
| `trajectoryDEG/` | `run_trajectory_gam_differential_gene_analysis` | Results CSV + heatmap/volcano/gene-curve PNGs |
| `sample_cluster/` | `cluster`, `proportion_test` | K-means cluster CSVs, embedding scatters, proportion test heatmaps |
| `raisin_results_*/` | `raisinfit` + `run_pairwise_tests` | Per-pair subdirs (`0_vs_1/…`) with volcano + results CSV, plus `summary_plots/` |
| `visualization/` | `visualization`, `visualize_multimodal_embedding` | Dendrogram, grouped scatter plots, multi-omics embedding panels |

## Two worked configs

=== "Minimal RNA from scratch"

    ```yaml
    output_dir: "/data/run/covid_rna"
    run_rna_pipeline: true
    use_gpu: true

    rna_count_data_path: "/data/test_RNA.h5ad"
    rna_sample_meta_path: "/data/sample_meta.csv"
    rna_sample_col: "sample"
    rna_celltype_col: "cell_type"
    rna_sample_level_batch_col: null

    rna_preprocessing: true
    rna_cell_type_cluster: true
    rna_derive_sample_embedding: true
    rna_sample_distance_calculation: true
    rna_trajectory_analysis: true
    rna_trajectory_dge: true
    rna_sample_cluster: true
    rna_visualize_data: true

    rna_min_cells: 500
    rna_min_genes: 500
    rna_pct_mito_cutoff: 20
    rna_num_cell_hvgs: 2000
    rna_leiden_cluster_resolution: 0.99
    rna_sample_hvg_number: 2000
    rna_sample_embedding_dimension: 10

    rna_trajectory_supervised: true
    rna_trajectory_col: "sev.level"
    rna_n_cca_pcs: 10
    rna_cca_pvalue: true

    rna_sample_distance_methods: ["cosine", "correlation"]
    rna_grouping_columns: ["sev.level"]
    rna_cluster_number: 4
    ```

=== "Resume from pseudobulk (RNA)"

    ```yaml
    output_dir: "/data/run/covid_rna"
    run_rna_pipeline: true
    use_gpu: true

    # skip preprocessing and embedding — already on disk
    rna_preprocessing: false
    rna_cell_type_cluster: false
    rna_derive_sample_embedding: false

    rna_adata_cell_path: "/data/run/covid_rna/rna/preprocess/adata_cell.h5ad"
    rna_adata_sample_path: "/data/run/covid_rna/rna/preprocess/adata_sample.h5ad"
    rna_pseudo_adata_path: "/data/run/covid_rna/rna/pseudobulk/pseudobulk_sample.h5ad"
    rna_sample_meta_path: "/data/sample_meta.csv"

    # run downstream only
    rna_sample_distance_calculation: true
    rna_trajectory_analysis: true
    rna_trajectory_dge: true
    rna_sample_cluster: true
    rna_trajectory_col: "sev.level"
    rna_sample_distance_methods: ["cosine"]
    rna_grouping_columns: ["sev.level"]
    rna_cluster_number: 4
    ```

## Pitfalls

!!! warning "Unknown keys abort the run"
    `SampleDisc.py` validates every config key against the wrapper signatures. A typo like `rna_leden_cluster_resolution` raises an immediate `TypeError` instead of being silently ignored.

!!! warning "Resume paths must match the upstream schema"
    If you disable `rna_preprocessing` you also must keep `rna_sample_col`, `rna_celltype_col`, and any `rna_*_batch_*` keys consistent with what was set when the files were originally written — the pseudobulk reads those columns from `.obs`.

!!! warning "use_gpu=true needs CUDA wheels"
    The `_linux` code path imports `rapids_singlecell`. If your environment only has CPU scanpy, set `use_gpu: false` and the wrapper will fall back to the CPU implementations.
