# Config key reference

<!-- GENERATED FILE: do not edit. Regenerate with docs/assets/gen_config_reference.py -->

Every key that `sampledisco -m complex --config` accepts, for sampledisco **0.10.0**: **311 keys** (Global 11, RNA 80, ATAC 85, Multi-omics 135).

Generated from the package by [`gen_config_reference.py`](../assets/gen_config_reference.py); do not edit by hand. **Default** is the `wrapper()` default, which a key left out of a Python call takes; the CLI requires every key except the 10 optional ones. **Demo** is the value in the template that `sampledisco --init-config` writes, in bold when it differs. **Consumed by** is the function and parameter the value reaches, or the stage it switches (`switch in X`). How the keys fit together is in the [Configuration guide](configuration.md).

!!! warning "15 keys are accepted but ignored"
    The CLI requires them, so leave them in the file, but changing them has no effect:

    - `rna_adata_sample_path`: `wrapper()` never reads it.
    - `rna_pseudo_adata_path`: `rna_wrapper()` accepts `sample_adata_path` but never reads it.
    - `atac_adata_sample_path`: `wrapper()` never reads it.
    - `atac_pseudo_adata_path`: `atac_wrapper()` accepts `sample_adata_path` but never reads it.
    - `multiomics_pseudobulk_h5ad_path`: `multiomics_wrapper()` accepts `sample_adata_path` but never reads it.
    - `multiomics_lsi_n_iter`: `multiomics_wrapper()` accepts `lsi_n_iter` but never reads it.
    - `multiomics_compression`: `multiomics_wrapper()` accepts `compression` but never reads it.
    - `multiomics_existing_cell_types`: `multiomics_wrapper()` accepts `existing_cell_types` but never reads it.
    - `multiomics_n_target_clusters`: `multiomics_wrapper()` accepts `n_target_clusters` but never reads it.
    - `multiomics_markers`: `multiomics_wrapper()` accepts `markers` but never reads it.
    - `multiomics_age_bin_size`: only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call.
    - `multiomics_age_column`: only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call.
    - `multiomics_plot_dendrogram_flag`: only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call.
    - `multiomics_plot_cell_type_proportions_pca_flag`: only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call.
    - `multiomics_plot_cell_type_expression_umap_flag`: only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call.

??? info "Template vs. signature check"
    - Keys in `wrapper()` but missing from the template (excluding optional keys): none.
    - Keys in the template that `wrapper()` does not accept: none.
    - Optional keys (the CLI fills in the default): `atac_batch_col`, `atac_save_cell_adata_after_embedding`, `atac_trajectory_test`, `multiomics_cell_level_batch_key`, `multiomics_sample_level_batch_col`, `multiomics_save_cell_adata_after_embedding`, `multiomics_trajectory_test`, `rna_batch_col`, `rna_save_cell_adata_after_embedding`, `rna_trajectory_test`.

## Global keys (11)

### Global: output

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `output_dir` | *(required)* | **`"./sampledisco_demo_output"`** | 20 functions: `preprocess`, `cell_types`, `run_autotune`, `compute_sample_embedding`, `sample_distance`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `run_pairwise_tests`, `visualization`, `visualize_multimodal_embedding`, `run_dimension_association_analysis`, `multiomics_preparation`, `cell_types_multiomics`, `propagate_cell_type`, `save_embedding_to_h5ad`, `warn_if_stale_embedding` |  |

### Global: pipeline selection

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `run_rna_pipeline` | `true` | `true` | switch in wrapper: rna_wrapper, downstream_analysis | Default `true`: set it to `false` explicitly for ATAC-only or multi-omics-only runs. |
| `run_atac_pipeline` | `false` | **`true`** | switch in wrapper: atac_wrapper, downstream_analysis |  |
| `run_multiomics_pipeline` | `false` | `false` | switch in wrapper: multiomics_wrapper, downstream_analysis | Template: *set true to run multi-omics; pre-wired to skip scGLUE training (see multiomics_integration below)* |

### Global: general settings

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `use_gpu` | `false` | `false` | `compute_sample_embedding:use_gpu`<br>switch in rna_wrapper: preprocess_gpu instead of preprocess<br>switch in rna_wrapper: cell_types_gpu instead of cell_types<br>switch in atac_wrapper: preprocess_gpu instead of preprocess<br>switch in atac_wrapper: cell_types_gpu instead of cell_types<br>`harmonize_xglue:use_gpu` | Off unless set: runs Harmony, neighbors, Leiden and UMAP on the GPU through RAPIDS (the embedding's k-means stays on scikit-learn). Falls back to CPU with a warning if RAPIDS cannot initialise; `gpu_available` in `sys_log/main_process_status.json` records the backend. Results differ from a CPU run; see [GPU acceleration](gpu.md). Template: *set true to use the RAPIDS GPU paths if installed; auto-falls back to CPU otherwise* |
| `initialization` | `true` | **`false`** | switch in wrapper | Default `true`. `true` deletes `{output_dir}/{rna,atac,multiomics}/result/` (a folder no stage writes) and rewrites `sys_log/main_process_status.json`; `false` reloads that JSON and prints it. Neither setting skips a stage: the stage flags alone decide what runs. |
| `verbose` | `true` | `true` | 12 functions: `preprocess`, `cell_types`, `run_autotune`, `compute_sample_embedding`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `raisinfit`, `visualization`, `visualize_multimodal_embedding`, `run_dimension_association_analysis` |  |
| `save_intermediate` | `true` | `true` | switch in multiomics_wrapper: save_embedding_to_h5ad | Only used when `multiomics_derive_sample_embedding`. Read only by the multi-omics pipeline (re-saves `preprocess/adata_sample.h5ad` after the Harmony pass and after the sample embedding). RNA and ATAC always write their files. |
| `large_data_need_extra_memory` | `false` | `false` | `reinitialize:managed_memory`<br>`reinitialize:pool_allocator` | Only used when `use_gpu`. GPU only: switches the RAPIDS memory pool to managed (oversubscribable) memory. No effect on CPU. |
| `continue_on_error` | `false` | `false` | switch in wrapper | `false`: the first failing modality aborts the run. `true`: the error is stored as `results['<modality>_error']`, later modalities still run, and the CLI still exits 1. |
| `random_state` | `42` | `42` | `compute_sample_embedding:seed`<br>`raisinfit:seed`<br>`run_pairwise_tests:random_state` | Seeds the RNA/ATAC sample embedding and RAISIN only. Preprocessing and cell typing use a fixed seed 42; multi-omics uses `multiomics_random_state`. |

## RNA keys (80)

### RNA: input and output paths

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_count_data_path` | `null` | **`"./data/test_RNA.h5ad"`** | `preprocess:h5ad_path` | Only used when `rna_preprocessing`. Raw counts in `.X`. Required whenever this pipeline runs, even when resuming. |
| `rna_output_dir` | `null` | `null` | 15 functions: `preprocess`, `cell_types`, `run_autotune`, `compute_sample_embedding`, `sample_distance`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `run_pairwise_tests`, `visualization`, `visualize_multimodal_embedding`, `run_dimension_association_analysis` | Template: *defaults to {output_dir}/rna* |

### RNA: pipeline control flags

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_preprocessing` | `true` | `true` | switch in rna_wrapper: preprocess |  |
| `rna_cell_type_cluster` | `true` | `true` | switch in rna_wrapper: cell_types |  |
| `rna_derive_sample_embedding` | `true` | `true` | switch in rna_wrapper: resolve_embedding_keys, run_autotune, compute_sample_embedding |  |
| `rna_sample_distance_calculation` | `true` | `true` | switch in downstream_analysis: sample_distance |  |
| `rna_trajectory_analysis` | `true` | `true` | switch in downstream_analysis: CCA_Call, TSCAN, cca_pvalue_test, run_trajectory_gam_differential_gene_analysis |  |
| `rna_trajectory_dge` | `true` | `true` | switch in downstream_analysis: run_trajectory_gam_differential_gene_analysis | Only used when `rna_trajectory_analysis`. |
| `rna_sample_cluster` | `true` | `true` | switch in downstream_analysis: cluster |  |
| `rna_proportion_test` | `false` | **`true`** | switch in downstream_analysis: proportion_test | Groups: `rna_cluster_differential_gene_group_col` if set, else the k-means clusters (needs `rna_sample_cluster: true`). |
| `rna_cluster_dge` | `false` | `false` | switch in downstream_analysis: raisinfit, run_pairwise_tests | RAISIN: memory-heavy (scales with cores and features). Off in the demo. Template: *RAISIN cluster-DGE is memory-heavy (multithreaded ComBat scales with #cores); off by default so the demo finishes on modest machines. Set true if you have ample RAM.* |
| `rna_visualize_data` | `true` | `true` | switch in downstream_analysis: visualization |  |
| `rna_dimension_association_analysis` | `false` | **`true`** | switch in downstream_analysis: run_dimension_association_analysis |  |

### RNA: input data paths

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_adata_cell_path` | `null` | `null` | `rna_wrapper:adata_path` | Resume input: read when `rna_preprocessing: false` (default `{output dir}/preprocess/adata_preprocessed.h5ad`). With `rna_derive_sample_embedding: false` the embedding is taken from this file's `uns['X_DR_sample']`. |
| `rna_adata_sample_path` | `null` | `null` | *nothing* | **Dead key:** `wrapper()` never reads it. |
| `rna_sample_meta_path` | `null` | `null` | `preprocess:sample_meta_path` | Only used when `rna_preprocessing`. Optional sample-level CSV, one row per sample, joined on the sample column. Its columns replace same-named `.obs` columns. |
| `rna_cell_meta_path` | `null` | `null` | `preprocess:cell_meta_path` | Only used when `rna_preprocessing`. Optional cell-level CSV with a `barcode` column matching `obs_names`. |
| `rna_pseudo_adata_path` | `null` | `null` | *nothing* | **Dead key:** `rna_wrapper()` accepts `sample_adata_path` but never reads it. |

### RNA: common column names

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_sample_col` | `"sample"` | `"sample"` | 9 functions: `build_sample_adata`, `preprocess`, `run_autotune`, `compute_sample_embedding`, `sample_distance`, `run_trajectory_gam_differential_gene_analysis`, `proportion_test`, `raisinfit`, `run_dimension_association_analysis` | A column literally named `sample` is safest. With a sample-metadata CSV, this column is kept and a `sample` copy is added. |
| `rna_batch_col` | `null` | **`"batch"`** | `resolve_batch_levels:batch_col` |  |
| `rna_sample_level_batch_col` | `null` | `null` | `resolve_batch_levels:sample_level` | A string or a list. Two or more columns give multi-covariate sample-level Harmony. Harmony runs only with at least 8 samples and more than one batch value. |
| `rna_celltype_col` | `"cell_type"` | `"cell_type"` | `cell_types:cell_type_column`<br>`run_autotune:celltype_col`<br>`compute_sample_embedding:celltype_col`<br>`sample_distance:cell_type_column`<br>`run_trajectory_gam_differential_gene_analysis:celltype_col`<br>`proportion_test:celltype_col` | Cell typing always writes its labels to `obs['cell_type']`, and once it has run the embedding and downstream steps use that column. This key is the label source with `rna_existing_cell_types: true`, and the column the embedding reads when `rna_cell_type_cluster: false`. Template: *existing labels to start from (existing_cell_types:true) and the label column used when cell typing is off. Cell typing always writes obs['cell_type'], which later stages then use. The demo .h5ad's precomputed labels are in 'celltype' — set this to "celltype" only if you flip existing_cell_types:true.* |

### RNA: preprocessing parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_min_cells` | `500` | **`100`** | `preprocess:min_cells` | Only used when `rna_preprocessing`. Used twice: genes detected in fewer cells are dropped, and samples with fewer cells are dropped. |
| `rna_min_genes` | `500` | **`200`** | `preprocess:min_genes` | Only used when `rna_preprocessing`. |
| `rna_pct_mito_cutoff` | `20` | `20.0` | `preprocess:pct_mito_cutoff` | Only used when `rna_preprocessing`. |
| `rna_exclude_genes` | `null` | `null` | `preprocess:exclude_genes` | Only used when `rna_preprocessing`. |
| `rna_num_cell_hvgs` | `2000` | `2000` | `preprocess:num_cell_hvgs` | Only used when `rna_preprocessing`. |
| `rna_cell_embedding_num_pcs` | `20` | `20` | `preprocess:cell_embedding_num_PCs`<br>`cell_types:cell_embedding_num_PCs` |  |
| `rna_num_harmony_iterations` | `30` | **`20`** | `preprocess:num_harmony_iterations` | Only used when `rna_preprocessing`. |
| `rna_cell_level_batch_key` | `null` | `null` | `resolve_batch_levels:cell_level` |  |

### RNA: cell type clustering parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_leiden_cluster_resolution` | `0.8` | `0.8` | `cell_types:leiden_cluster_resolution` | Only used when `rna_cell_type_cluster`. |
| `rna_cell_embedding_column` | `null` | `null` | `cell_types:cell_embedding_column`<br>`run_autotune:comp_emb_key`<br>`run_autotune:rmd_emb_key`<br>`compute_sample_embedding:comp_emb_key`<br>`compute_sample_embedding:rmd_emb_key` |  |
| `rna_existing_cell_types` | `false` | `false` | `cell_types:existing_cell_types` | Only used when `rna_cell_type_cluster`. Use labels already in `obs[rna_celltype_col]`; they are copied to `obs['cell_type']` (and merged down to `rna_n_target_cell_clusters` if that is smaller, with a warning). |
| `rna_n_target_cell_clusters` | `null` | **`8`** | `cell_types:n_target_clusters` | Only used when `rna_cell_type_cluster`. |
| `rna_umap` | `false` | **`true`** | `cell_types:umap` | Only used when `rna_cell_type_cluster`. |

### RNA: sample embedding parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_sample_embedding_medium_K` | `120` | `120` | `run_autotune:medium_K`<br>`compute_sample_embedding:medium_K` | Only used when `rna_derive_sample_embedding`. |
| `rna_sample_embedding_fine_K` | `300` | `300` | `run_autotune:fine_K`<br>`compute_sample_embedding:fine_K` | Only used when `rna_derive_sample_embedding`. |
| `rna_sample_embedding_rmd_dim` | `8` | `8` | `run_autotune:rmd_dim`<br>`compute_sample_embedding:rmd_dim_per_cluster` | Only used when `rna_derive_sample_embedding`. |
| `rna_sample_embedding_use_clr` | `false` | `false` | `compute_sample_embedding:use_clr` | Only used when `not rna_autotune_enable` and `rna_derive_sample_embedding`. Ignored when `rna_autotune_enable: true`. |
| `rna_sample_embedding_use_rmd` | `true` | `true` | `compute_sample_embedding:use_rmd` | Only used when `not rna_autotune_enable` and `rna_derive_sample_embedding`. Ignored when `rna_autotune_enable: true`. |
| `rna_sample_embedding_block_weights` | `null` | `null` | `compute_sample_embedding:block_weights` | Only used when `not rna_autotune_enable` and `rna_derive_sample_embedding`. Ignored when `rna_autotune_enable: true`. |
| `rna_sample_embedding_rmd_weight` | `"equal"` | `"equal"` | `compute_sample_embedding:rmd_weight` | Only used when `not rna_autotune_enable` and `rna_derive_sample_embedding`. Validated by the CLI: `"equal"` or a positive number. Ignored when `rna_autotune_enable: true`. Template: *RMD block weight α: "equal" = same energy as the three composition blocks together (α² = w_A1² + w_A2² + w_A3²); a number fixes α (0.60 = previous default)* |
| `rna_sample_embedding_pca_components` | `10` | `10` | `run_autotune:pca_components`<br>`compute_sample_embedding:pca_components` | Only used when `rna_derive_sample_embedding`. Capped at n_units - 1 (8 demo samples give 7 PCs). |
| `rna_sample_embedding_batch_method` | `"harmony"` | `"harmony"` | `run_autotune:batch_method`<br>`compute_sample_embedding:batch_method` | Only used when `rna_derive_sample_embedding`. `"harmony"`, `"linear"` or `"none"`; any other value is rejected by the CLI and by `compute_sample_embedding`. |

### RNA: autotune parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_autotune_enable` | `false` | `false` | switch in rna_wrapper: run_autotune | Only used when `rna_derive_sample_embedding`. Replaces `compute_sample_embedding` with `run_autotune` (CPU only). Without a batch column and without `rna_autotune_grouping_col` no search runs and alpha stays "equal". |
| `rna_autotune_search` | `"bayesian"` | `"bayesian"` | `run_autotune:search` | Only used when `rna_autotune_enable` and `rna_derive_sample_embedding`. |
| `rna_autotune_scoring` | `"auto"` | `"auto"` | `run_autotune:scoring` | Only used when `rna_autotune_enable` and `rna_derive_sample_embedding`. |
| `rna_autotune_scope` | `"alpha_only"` | `"alpha_only"` | `run_autotune:scope` | Only used when `rna_autotune_enable` and `rna_derive_sample_embedding`. |
| `rna_autotune_alpha_bounds` | `[0.1, 100.0]` | `[0.1, 100.0]` | `run_autotune:alpha_bounds` | Only used when `rna_autotune_enable` and `rna_derive_sample_embedding`. |
| `rna_save_cell_adata_after_embedding` | `true` | `true` | `run_autotune:save_cell_adata`<br>`compute_sample_embedding:save_cell_adata` | Only used when `rna_derive_sample_embedding`. Optional: may be left out of a config (default `true`). Template: *store .uns['X_DR_sample'] in the cell-level h5ad (updated in place; false = keep only the CSV)* |
| `rna_autotune_grouping_col` | `null` | `null` | `run_autotune:grouping_col` | Only used when `rna_autotune_enable` and `rna_derive_sample_embedding`. |

### RNA: trajectory analysis parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_n_cca_pcs` | `2` | **`5`** | `CCA_Call:n_components`<br>`cca_pvalue_test:n_pcs` | Only used when `rna_trajectory_analysis` and `rna_trajectory_supervised`. Capped at the number of embedding PCs. Keep it well below the sample count: with PCs close to n_samples - 1 CCA fits any labelling perfectly. |
| `rna_trajectory_col` | `"sev.level"` | `"sev.level"` | `CCA_Call:trajectory_col`<br>`cca_pvalue_test:trajectory_col`<br>`run_trajectory_gam_differential_gene_analysis:anchor_col` | Only used when `rna_trajectory_analysis`. **Change for your data.** Default `sev.level`. Coerced to numeric for CCA: text labels become NaN and those samples are dropped. |
| `rna_trajectory_supervised` | `false` | **`true`** | `run_trajectory_gam_differential_gene_analysis:anchor_col`<br>switch in downstream_analysis: CCA_Call, cca_pvalue_test | Only used when `rna_trajectory_analysis`. |
| `rna_trajectory_visualization_label` | `null` | **`["sev.level"]`** | `TSCAN:grouping_columns` | Only used when `not rna_trajectory_supervised` and `rna_trajectory_analysis`. **Change for your data.** `null` becomes `["sev.level"]` (TSCAN plot colouring). |
| `rna_cca_pvalue` | `false` | `false` | switch in downstream_analysis: cca_pvalue_test | Only used when `rna_trajectory_analysis` and `rna_trajectory_supervised`. Permutation test with the same number of PCs as the observed score. When the PCs approach n_samples - 1 every permuted score is ~1: the step warns and the p-value is meaningless. |
| `rna_tscan_origin` | `null` | `null` | `TSCAN:origin` | Only used when `not rna_trajectory_supervised` and `rna_trajectory_analysis`. An integer: the index of the endpoint cluster to start from. |
| `rna_tscan_n_clusters` | `null` | **`3`** | `TSCAN:n_clusters` | Only used when `not rna_trajectory_supervised` and `rna_trajectory_analysis`. Template: *int; null lets TSCAN auto-select via BIC (2-20)* |
| `rna_tscan_pseudotime_mode` | `"rank"` | `"rank"` | `TSCAN:pseudotime_mode` | Only used when `not rna_trajectory_supervised` and `rna_trajectory_analysis`. `"rank"` or `"distance"`; anything else raises. Template: *"rank" or "distance"* |

### RNA: sample distance parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_sample_distance_methods` | `null` | **`["cosine", "correlation"]`** | `sample_distance:method` | Only used when `rna_sample_distance_calculation`. SciPy metrics on the embedding (`euclidean`, `cosine`, `correlation`, … 12 in all) and `EMD`, `chi_square`, `jensen_shannon` on the cell-type composition; unknown names are skipped with a warning. |
| `rna_grouping_columns` | `null` | **`["sev.level"]`** | `sample_distance:grouping_columns`<br>`visualization:grouping_columns` | **Change for your data.** `null` becomes `["sev.level"]`. Every column must exist in the sample metadata, or the visualization step raises `ValueError: Missing required columns`. |
| `rna_summary_sample_csv_path` | `null` | `null` | `sample_distance:summary_csv_path` | Only used when `rna_sample_distance_calculation`. Written by the distance step and deleted after the trajectory step: always for the default path, and for a path you supply only if the run created it. |

### RNA: trajectory DGE parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_fdr_threshold` | `0.05` | `0.05` | `run_trajectory_gam_differential_gene_analysis:fdr_threshold` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. |
| `rna_effect_size_threshold` | `1` | `1.0` | `run_trajectory_gam_differential_gene_analysis:effect_size_threshold` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. |
| `rna_top_n_genes` | `100` | `100` | `run_trajectory_gam_differential_gene_analysis:top_n_genes` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. |
| `rna_trajectory_diff_gene_covariate` | `null` | `null` | `run_trajectory_gam_differential_gene_analysis:covariate_columns` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. |
| `rna_trajectory_test` | `"gam"` | `"gam"` | `run_trajectory_gam_differential_gene_analysis:test` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. Template: *"gam" or "linear"; linear suits small cohorts* |
| `rna_num_splines` | `5` | `5` | `run_trajectory_gam_differential_gene_analysis:num_splines` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. |
| `rna_spline_order` | `3` | **`2`** | `run_trajectory_gam_differential_gene_analysis:spline_order` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. A gene is fitted only if it has a value in at least order + splines + 2 samples (after automatic reduction). With 8 samples use 2; order 3 fits nothing. |
| `rna_visualization_gene_list` | `null` | `null` | `run_trajectory_gam_differential_gene_analysis:visualization_gene_list` | Only used when `rna_trajectory_analysis` and `rna_trajectory_dge`. |

### RNA: sample clustering parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_cluster_number` | `4` | **`3`** | `cluster:number_of_clusters` | Only used when `rna_sample_cluster`. |
| `rna_cluster_differential_gene_group_col` | `null` | **`"sev.level"`** | `proportion_test:group_col`<br>`raisinfit:group_col`<br>switch in downstream_analysis: proportion_test<br>switch in downstream_analysis: raisinfit, run_pairwise_tests | When set, the proportion test and RAISIN compare these groups **instead of** the k-means sample clusters (the demo sets `sev.level`). |

### RNA: visualization parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_age_bin_size` | `null` | `null` | `visualization:age_bin_size` | Only used when `rna_visualize_data`. |
| `rna_age_column` | `"age"` | `"age"` | `visualization:age_column` | Only used when `rna_visualize_data`. |
| `rna_plot_dendrogram_flag` | `true` | `true` | `visualization:plot_dendrogram_flag` | Only used when `rna_visualize_data`. |
| `rna_plot_cell_type_proportions_pca_flag` | `false` | **`true`** | `visualization:plot_cell_type_proportions_pca_flag` | Only used when `rna_visualize_data`. |
| `rna_plot_cell_type_expression_umap_flag` | `false` | **`true`** | `visualization:plot_cell_type_expression_umap_flag` | Only used when `rna_visualize_data`. |

### RNA: dimension association analysis parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `rna_association_continuous_cols` | `null` | `null` | `run_dimension_association_analysis:continuous_cols` | Only used when `rna_dimension_association_analysis`. |
| `rna_association_categorical_cols` | `null` | `null` | `run_dimension_association_analysis:categorical_cols` | Only used when `rna_dimension_association_analysis`. |
| `rna_association_n_permutations` | `999` | **`99`** | `run_dimension_association_analysis:n_permutations` | Only used when `rna_dimension_association_analysis`. |

## ATAC keys (85)

### ATAC: input and output paths

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_count_data_path` | `null` | **`"./data/test_ATAC.h5ad"`** | `preprocess:h5ad_path` | Only used when `atac_preprocessing`. Raw counts in `.X`. Required whenever this pipeline runs, even when resuming. |
| `atac_output_dir` | `null` | `null` | 15 functions: `preprocess`, `cell_types`, `run_autotune`, `compute_sample_embedding`, `sample_distance`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `run_pairwise_tests`, `visualization`, `visualize_multimodal_embedding`, `run_dimension_association_analysis` | Template: *defaults to {output_dir}/atac* |

### ATAC: pipeline control flags

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_preprocessing` | `true` | `true` | switch in atac_wrapper: preprocess |  |
| `atac_cell_type_cluster` | `true` | `true` | switch in atac_wrapper: cell_types |  |
| `atac_derive_sample_embedding` | `true` | `true` | switch in atac_wrapper: resolve_embedding_keys, run_autotune, compute_sample_embedding |  |
| `atac_sample_distance_calculation` | `true` | `true` | switch in downstream_analysis: sample_distance |  |
| `atac_trajectory_analysis` | `true` | `true` | switch in downstream_analysis: CCA_Call, TSCAN, cca_pvalue_test, run_trajectory_gam_differential_gene_analysis |  |
| `atac_trajectory_dge` | `true` | **`false`** | switch in downstream_analysis: run_trajectory_gam_differential_gene_analysis | Only used when `atac_trajectory_analysis`. Template: *8 samples × 4 splines: GAM cannot fit (too_few_samples for every gene)* |
| `atac_sample_cluster` | `true` | `true` | switch in downstream_analysis: cluster |  |
| `atac_proportion_test` | `false` | **`true`** | switch in downstream_analysis: proportion_test | Groups: `atac_cluster_differential_gene_group_col` if set, else the k-means clusters (needs `atac_sample_cluster: true`). |
| `atac_cluster_dge` | `false` | `false` | switch in downstream_analysis: raisinfit, run_pairwise_tests | RAISIN: memory-heavy (scales with cores and features). Off in the demo. Template: *disabled: 29K cells × 230K peaks × 48 cores RAISIN OOM-kills on CPU* |
| `atac_visualize_data` | `true` | `true` | switch in downstream_analysis: visualization |  |
| `atac_dimension_association_analysis` | `false` | **`true`** | switch in downstream_analysis: run_dimension_association_analysis |  |

### ATAC: input data paths

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_adata_cell_path` | `null` | `null` | `atac_wrapper:adata_path` | Resume input: read when `atac_preprocessing: false` (default `{output dir}/preprocess/adata_preprocessed.h5ad`). With `atac_derive_sample_embedding: false` the embedding is taken from this file's `uns['X_DR_sample']`. |
| `atac_adata_sample_path` | `null` | `null` | *nothing* | **Dead key:** `wrapper()` never reads it. |
| `atac_sample_meta_path` | `null` | `null` | `preprocess:sample_meta_path` | Only used when `atac_preprocessing`. Optional sample-level CSV, one row per sample, joined on the sample column. Its columns replace same-named `.obs` columns. |
| `atac_cell_meta_path` | `null` | `null` | `preprocess:cell_meta_path` | Only used when `atac_preprocessing`. Optional cell-level CSV with a `barcode` column matching `obs_names`. |
| `atac_pseudo_adata_path` | `null` | `null` | *nothing* | **Dead key:** `atac_wrapper()` accepts `sample_adata_path` but never reads it. |

### ATAC: common column names

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_sample_col` | `"sample"` | `"sample"` | 9 functions: `build_sample_adata`, `preprocess`, `run_autotune`, `compute_sample_embedding`, `sample_distance`, `run_trajectory_gam_differential_gene_analysis`, `proportion_test`, `raisinfit`, `run_dimension_association_analysis` | A column literally named `sample` is safest. With a sample-metadata CSV, this column is kept and a `sample` copy is added. |

### ATAC: batch levels as for RNA

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_batch_col` | `null` | `null` | `resolve_batch_levels:batch_col` | Template: *default for both batch levels, as for RNA* |
| `atac_sample_level_batch_col` | `null` | `null` | `resolve_batch_levels:sample_level` | A string or a list. Two or more columns give multi-covariate sample-level Harmony. Harmony runs only with at least 8 samples and more than one batch value. |
| `atac_celltype_col` | `"cell_type"` | `"cell_type"` | `cell_types:cell_type_column`<br>`run_autotune:celltype_col`<br>`compute_sample_embedding:celltype_col`<br>`sample_distance:cell_type_column`<br>`run_trajectory_gam_differential_gene_analysis:celltype_col`<br>`proportion_test:celltype_col` | Cell typing always writes its labels to `obs['cell_type']`, and once it has run the embedding and downstream steps use that column. This key is the label source with `atac_existing_cell_types: true`, and the column the embedding reads when `atac_cell_type_cluster: false`. |
| `atac_cell_embedding_column` | `null` | `null` | `cell_types:cell_embedding_column`<br>`run_autotune:comp_emb_key`<br>`run_autotune:rmd_emb_key`<br>`compute_sample_embedding:comp_emb_key`<br>`compute_sample_embedding:rmd_emb_key` |  |

### ATAC: ATAC-specific preprocessing parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_min_cells` | `1` | **`10`** | `preprocess:min_cells` | Only used when `atac_preprocessing`. |
| `atac_min_features` | `2000` | **`200`** | `preprocess:min_features` | Only used when `atac_preprocessing`. |
| `atac_max_features` | `15000` | **`50000`** | `preprocess:max_features` | Only used when `atac_preprocessing`. |
| `atac_min_cells_per_sample` | `1` | **`50`** | `preprocess:min_cells_per_sample` | Only used when `atac_preprocessing`. |
| `atac_exclude_features` | `null` | `null` | `preprocess:exclude_features` | Only used when `atac_preprocessing`. |
| `atac_cell_level_batch_key` | `null` | `null` | `resolve_batch_levels:cell_level` |  |
| `atac_doublet_detection` | `true` | **`false`** | `preprocess:doublet_detection` | Only used when `atac_preprocessing`. |
| `atac_num_cell_hvfs` | `50000` | **`10000`** | `preprocess:num_cell_hvfs` | Only used when `atac_preprocessing`. |
| `atac_cell_embedding_num_pcs` | `50` | **`30`** | `preprocess:cell_embedding_num_PCs`<br>`cell_types:cell_embedding_num_PCs` |  |
| `atac_num_harmony_iterations` | `30` | **`20`** | `preprocess:num_harmony_iterations` | Only used when `atac_preprocessing`. |
| `atac_tfidf_scale_factor` | `10000.0` | `10000.0` | `preprocess:tfidf_scale_factor` | Only used when `atac_preprocessing`. |
| `atac_log_transform` | `true` | `true` | `preprocess:log_transform` | Only used when `atac_preprocessing`. |
| `atac_drop_first_lsi` | `true` | `true` | `preprocess:drop_first_lsi` | Only used when `atac_preprocessing`. |

### ATAC: cell type clustering parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_leiden_cluster_resolution` | `0.8` | `0.8` | `cell_types:leiden_cluster_resolution` | Only used when `atac_cell_type_cluster`. |
| `atac_existing_cell_types` | `false` | `false` | `cell_types:existing_cell_types` | Only used when `atac_cell_type_cluster`. Use labels already in `obs[atac_celltype_col]`; they are copied to `obs['cell_type']` (and merged down to `atac_n_target_cell_clusters` if that is smaller, with a warning). |
| `atac_n_target_cell_clusters` | `null` | **`8`** | `cell_types:n_target_clusters` | Only used when `atac_cell_type_cluster`. |
| `atac_umap` | `false` | **`true`** | `cell_types:umap` | Only used when `atac_cell_type_cluster`. |

### ATAC: sample embedding parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_sample_embedding_medium_K` | `120` | `120` | `run_autotune:medium_K`<br>`compute_sample_embedding:medium_K` | Only used when `atac_derive_sample_embedding`. |
| `atac_sample_embedding_fine_K` | `300` | `300` | `run_autotune:fine_K`<br>`compute_sample_embedding:fine_K` | Only used when `atac_derive_sample_embedding`. |
| `atac_sample_embedding_rmd_dim` | `8` | `8` | `run_autotune:rmd_dim`<br>`compute_sample_embedding:rmd_dim_per_cluster` | Only used when `atac_derive_sample_embedding`. |
| `atac_sample_embedding_use_clr` | `false` | `false` | `compute_sample_embedding:use_clr` | Only used when `atac_derive_sample_embedding` and `not atac_autotune_enable`. Ignored when `atac_autotune_enable: true`. |
| `atac_sample_embedding_use_rmd` | `true` | `true` | `compute_sample_embedding:use_rmd` | Only used when `atac_derive_sample_embedding` and `not atac_autotune_enable`. Ignored when `atac_autotune_enable: true`. |
| `atac_sample_embedding_block_weights` | `null` | `null` | `compute_sample_embedding:block_weights` | Only used when `atac_derive_sample_embedding` and `not atac_autotune_enable`. Ignored when `atac_autotune_enable: true`. |
| `atac_sample_embedding_rmd_weight` | `"equal"` | `"equal"` | `compute_sample_embedding:rmd_weight` | Only used when `atac_derive_sample_embedding` and `not atac_autotune_enable`. Validated by the CLI: `"equal"` or a positive number. Ignored when `atac_autotune_enable: true`. Template: *RMD block weight α: "equal" = same energy as the three composition blocks together (α² = w_A1² + w_A2² + w_A3²); a number fixes α (0.60 = previous default)* |
| `atac_sample_embedding_pca_components` | `10` | `10` | `run_autotune:pca_components`<br>`compute_sample_embedding:pca_components` | Only used when `atac_derive_sample_embedding`. Capped at n_units - 1 (8 demo samples give 7 PCs). |
| `atac_sample_embedding_batch_method` | `"harmony"` | `"harmony"` | `run_autotune:batch_method`<br>`compute_sample_embedding:batch_method` | Only used when `atac_derive_sample_embedding`. `"harmony"`, `"linear"` or `"none"`; any other value is rejected by the CLI and by `compute_sample_embedding`. |

### ATAC: autotune parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_autotune_enable` | `false` | `false` | switch in atac_wrapper: run_autotune | Only used when `atac_derive_sample_embedding`. Replaces `compute_sample_embedding` with `run_autotune` (CPU only). Without a batch column and without `atac_autotune_grouping_col` no search runs and alpha stays "equal". |
| `atac_autotune_search` | `"bayesian"` | `"bayesian"` | `run_autotune:search` | Only used when `atac_autotune_enable` and `atac_derive_sample_embedding`. |
| `atac_autotune_scoring` | `"auto"` | `"auto"` | `run_autotune:scoring` | Only used when `atac_autotune_enable` and `atac_derive_sample_embedding`. |
| `atac_autotune_scope` | `"alpha_only"` | `"alpha_only"` | `run_autotune:scope` | Only used when `atac_autotune_enable` and `atac_derive_sample_embedding`. |
| `atac_autotune_alpha_bounds` | `[0.1, 100.0]` | `[0.1, 100.0]` | `run_autotune:alpha_bounds` | Only used when `atac_autotune_enable` and `atac_derive_sample_embedding`. |
| `atac_save_cell_adata_after_embedding` | `true` | `true` | `run_autotune:save_cell_adata`<br>`compute_sample_embedding:save_cell_adata` | Only used when `atac_derive_sample_embedding`. Optional: may be left out of a config (default `true`). Template: *store .uns['X_DR_sample'] in the cell-level h5ad (updated in place; false = keep only the CSV)* |
| `atac_autotune_grouping_col` | `null` | `null` | `run_autotune:grouping_col` | Only used when `atac_autotune_enable` and `atac_derive_sample_embedding`. |

### ATAC: trajectory analysis parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_n_cca_pcs` | `2` | **`5`** | `CCA_Call:n_components`<br>`cca_pvalue_test:n_pcs` | Only used when `atac_trajectory_analysis` and `atac_trajectory_supervised`. Capped at the number of embedding PCs. Keep it well below the sample count: with PCs close to n_samples - 1 CCA fits any labelling perfectly. |
| `atac_trajectory_col` | `"sev.level"` | `"sev.level"` | `CCA_Call:trajectory_col`<br>`cca_pvalue_test:trajectory_col`<br>`run_trajectory_gam_differential_gene_analysis:anchor_col` | Only used when `atac_trajectory_analysis`. **Change for your data.** Default `sev.level`. Coerced to numeric for CCA: text labels become NaN and those samples are dropped. |
| `atac_trajectory_supervised` | `true` | `true` | `run_trajectory_gam_differential_gene_analysis:anchor_col`<br>switch in downstream_analysis: CCA_Call, cca_pvalue_test | Only used when `atac_trajectory_analysis`. |
| `atac_trajectory_visualization_label` | `null` | **`["sev.level"]`** | `TSCAN:grouping_columns` | Only used when `atac_trajectory_analysis` and `not atac_trajectory_supervised`. **Change for your data.** `null` becomes `["sev.level"]` (TSCAN plot colouring). |
| `atac_cca_pvalue` | `false` | `false` | switch in downstream_analysis: cca_pvalue_test | Only used when `atac_trajectory_analysis` and `atac_trajectory_supervised`. Permutation test with the same number of PCs as the observed score. When the PCs approach n_samples - 1 every permuted score is ~1: the step warns and the p-value is meaningless. |
| `atac_tscan_origin` | `null` | `null` | `TSCAN:origin` | Only used when `atac_trajectory_analysis` and `not atac_trajectory_supervised`. An integer: the index of the endpoint cluster to start from. |
| `atac_tscan_n_clusters` | `null` | **`3`** | `TSCAN:n_clusters` | Only used when `atac_trajectory_analysis` and `not atac_trajectory_supervised`. Template: *int; null lets TSCAN auto-select via BIC (2-20)* |
| `atac_tscan_pseudotime_mode` | `"rank"` | `"rank"` | `TSCAN:pseudotime_mode` | Only used when `atac_trajectory_analysis` and `not atac_trajectory_supervised`. `"rank"` or `"distance"`; anything else raises. Template: *"rank" or "distance"* |

### ATAC: sample distance parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_sample_distance_methods` | `null` | **`["cosine", "correlation"]`** | `sample_distance:method` | Only used when `atac_sample_distance_calculation`. SciPy metrics on the embedding (`euclidean`, `cosine`, `correlation`, … 12 in all) and `EMD`, `chi_square`, `jensen_shannon` on the cell-type composition; unknown names are skipped with a warning. |
| `atac_grouping_columns` | `null` | **`["sev.level"]`** | `sample_distance:grouping_columns`<br>`visualization:grouping_columns` | **Change for your data.** `null` becomes `["sev.level"]`. Every column must exist in the sample metadata, or the visualization step raises `ValueError: Missing required columns`. |
| `atac_summary_sample_csv_path` | `null` | `null` | `sample_distance:summary_csv_path` | Only used when `atac_sample_distance_calculation`. Written by the distance step and deleted after the trajectory step: always for the default path, and for a path you supply only if the run created it. |

### ATAC: trajectory DGE parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_fdr_threshold` | `0.05` | `0.05` | `run_trajectory_gam_differential_gene_analysis:fdr_threshold` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. |
| `atac_effect_size_threshold` | `1.0` | `1.0` | `run_trajectory_gam_differential_gene_analysis:effect_size_threshold` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. |
| `atac_top_n_genes` | `100` | `100` | `run_trajectory_gam_differential_gene_analysis:top_n_genes` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. |
| `atac_trajectory_diff_gene_covariate` | `null` | `null` | `run_trajectory_gam_differential_gene_analysis:covariate_columns` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. |
| `atac_trajectory_test` | `"gam"` | `"gam"` | `run_trajectory_gam_differential_gene_analysis:test` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. Template: *"gam" or "linear"; linear suits small cohorts* |
| `atac_num_splines` | `5` | `5` | `run_trajectory_gam_differential_gene_analysis:num_splines` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. |
| `atac_spline_order` | `3` | `3` | `run_trajectory_gam_differential_gene_analysis:spline_order` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. A gene is fitted only if it has a value in at least order + splines + 2 samples (after automatic reduction). With 8 samples use 2; order 3 fits nothing. |
| `atac_visualization_gene_list` | `null` | `null` | `run_trajectory_gam_differential_gene_analysis:visualization_gene_list` | Only used when `atac_trajectory_analysis` and `atac_trajectory_dge`. |

### ATAC: sample clustering parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_cluster_number` | `4` | **`3`** | `cluster:number_of_clusters` | Only used when `atac_sample_cluster`. |
| `atac_cluster_differential_gene_group_col` | `null` | **`"sev.level"`** | `proportion_test:group_col`<br>`raisinfit:group_col`<br>switch in downstream_analysis: proportion_test<br>switch in downstream_analysis: raisinfit, run_pairwise_tests | When set, the proportion test and RAISIN compare these groups **instead of** the k-means sample clusters (the demo sets `sev.level`). |

### ATAC: visualization parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_age_bin_size` | `null` | `null` | `visualization:age_bin_size` | Only used when `atac_visualize_data`. |
| `atac_age_column` | `"age"` | `"age"` | `visualization:age_column` | Only used when `atac_visualize_data`. |
| `atac_plot_dendrogram_flag` | `true` | `true` | `visualization:plot_dendrogram_flag` | Only used when `atac_visualize_data`. |
| `atac_plot_cell_type_proportions_pca_flag` | `false` | `false` | `visualization:plot_cell_type_proportions_pca_flag` | Only used when `atac_visualize_data`. |
| `atac_plot_cell_type_expression_umap_flag` | `false` | `false` | `visualization:plot_cell_type_expression_umap_flag` | Only used when `atac_visualize_data`. |

### ATAC: dimension association analysis parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `atac_association_continuous_cols` | `null` | `null` | `run_dimension_association_analysis:continuous_cols` | Only used when `atac_dimension_association_analysis`. |
| `atac_association_categorical_cols` | `null` | `null` | `run_dimension_association_analysis:categorical_cols` | Only used when `atac_dimension_association_analysis`. |
| `atac_association_n_permutations` | `999` | **`99`** | `run_dimension_association_analysis:n_permutations` | Only used when `atac_dimension_association_analysis`. |

## Multi-omics keys (135)

### Multi-omics: input and output paths

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_rna_file` | `null` | **`"./data/test_RNA.h5ad"`** | `multiomics_preparation:rna_file` | Only used when `multiomics_integration`. Raw counts for scGLUE training: needed only with `multiomics_integration: true`. |
| `multiomics_atac_file` | `null` | **`"./data/test_ATAC.h5ad"`** | `multiomics_preparation:atac_file` | Only used when `multiomics_integration`. Raw counts for scGLUE training: needed only with `multiomics_integration: true`. |
| `multiomics_output_dir` | `null` | `null` | 17 functions: `multiomics_preparation`, `cell_types_multiomics`, `propagate_cell_type`, `run_autotune`, `compute_sample_embedding`, `save_embedding_to_h5ad`, `warn_if_stale_embedding`, `sample_distance`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `run_pairwise_tests`, `visualize_multimodal_embedding`, `run_dimension_association_analysis` | Template: *defaults to {output_dir}/multiomics* |

### Multi-omics: pipeline control flags

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_integration` | `true` | **`false`** | switch in multiomics_wrapper: multiomics_preparation | Template: *false -> SKIP scGLUE training, load the pre-integrated file below. Set true to train from scratch (needs bedtools) + set the resume path to null.* |
| `multiomics_sample_distance_calculation` | `true` | `true` | switch in downstream_analysis: sample_distance |  |
| `multiomics_trajectory_analysis` | `true` | `true` | switch in downstream_analysis: CCA_Call, TSCAN, cca_pvalue_test, run_trajectory_gam_differential_gene_analysis |  |
| `multiomics_trajectory_dge` | `true` | **`false`** | switch in downstream_analysis: run_trajectory_gam_differential_gene_analysis | Only used when `multiomics_trajectory_analysis`. Template: *16 (sample, modality) units: GAM cannot fit on small test* |
| `multiomics_sample_cluster` | `true` | `true` | switch in downstream_analysis: cluster |  |
| `multiomics_proportion_test` | `false` | **`true`** | switch in downstream_analysis: proportion_test | Groups: `multiomics_cluster_differential_gene_group_col` if set, else the k-means clusters (needs `multiomics_sample_cluster: true`). |
| `multiomics_cluster_dge` | `false` | `false` | switch in downstream_analysis: raisinfit, run_pairwise_tests | RAISIN: memory-heavy (scales with cores and features). Off in the demo. Template: *disabled: same RAISIN OOM risk as ATAC (large peak matrix)* |
| `multiomics_visualize_embedding` | `true` | `true` | switch in downstream_analysis: visualize_multimodal_embedding |  |
| `multiomics_dimension_association_analysis` | `false` | **`true`** | switch in downstream_analysis: run_dimension_association_analysis |  |

### Multi-omics: GLUE sub-pipeline flags

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_run_glue_preprocessing` | `true` | `true` | `multiomics_preparation:run_preprocessing` | Only used when `multiomics_integration`. |
| `multiomics_run_glue_training` | `true` | `true` | `multiomics_preparation:run_training` | Only used when `multiomics_integration`. |
| `multiomics_run_glue_merge` | `true` | `true` | `multiomics_preparation:run_merge` | Only used when `multiomics_integration`. Template: *build embedding-only union (adata_sample.h5ad)* |
| `multiomics_run_glue_preprocess_per_modality` | `true` | `true` | `multiomics_preparation:run_preprocess_per_modality` | Only used when `multiomics_integration`. Template: *per-modality QC h5ads for DGE/RAISIN* |
| `multiomics_cell_type_cluster` | `true` | `true` | switch in multiomics_wrapper: cell_types_multiomics, propagate_cell_type | Labels RNA cells by Leiden and transfers them to ATAC. Modality values must be `RNA` / `ATAC`. |
| `multiomics_run_glue_visualization` | `true` | `true` | `multiomics_preparation:run_visualization`<br>`cell_types_multiomics:generate_plots` |  |

### Multi-omics: input data paths

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_integrated_h5ad_path` | `null` | **`"./data/test_multiomics_integrated.h5ad"`** | switch in multiomics_wrapper<br>switch in wrapper | Input only: read when `{multiomics output}/preprocess/adata_sample.h5ad` does not exist yet; every update goes to that output file. Template: *scGLUE-integrated object from Zenodo; used when multiomics_integration: false (skips training)* |
| `multiomics_pseudobulk_h5ad_path` | `null` | `null` | *nothing* | **Dead key:** `multiomics_wrapper()` accepts `sample_adata_path` but never reads it. |
| `multiomics_rna_sample_meta_file` | `null` | `null` | `multiomics_preparation:rna_sample_meta_file` | Only used when `multiomics_integration`. Optional sample-level CSV, one row per sample, joined on the sample column. Its columns replace same-named `.obs` columns. |
| `multiomics_atac_sample_meta_file` | `null` | `null` | `multiomics_preparation:atac_sample_meta_file` | Only used when `multiomics_integration`. Optional sample-level CSV, one row per sample, joined on the sample column. Its columns replace same-named `.obs` columns. |
| `multiomics_additional_hvg_file` | `null` | `null` | `multiomics_preparation:additional_hvg_file` | Only used when `multiomics_integration`. |

### Multi-omics: common column names

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_rna_sample_column` | `"sample"` | `"sample"` | `multiomics_preparation:rna_sample_column` | Only used when `multiomics_integration`. |
| `multiomics_atac_sample_column` | `"sample"` | `"sample"` | `multiomics_preparation:atac_sample_column` | Only used when `multiomics_integration`. |
| `multiomics_sample_col` | `"sample"` | `"sample"` | 10 functions: `build_sample_adata`, `multiomics_preparation`, `harmonize_xglue`, `run_autotune`, `compute_sample_embedding`, `sample_distance`, `run_trajectory_gam_differential_gene_analysis`, `proportion_test`, `raisinfit`, `run_dimension_association_analysis` | A column literally named `sample` is safest. With a sample-metadata CSV, this column is kept and a `sample` copy is added. |
| `multiomics_batch_col` | `null` | `null` | `resolve_batch_levels:batch_col` | A single column name. It is scGLUE's `batch_key`, the `harmonize_xglue` batch, and the sample-level Harmony batch. Template: *default for both batch levels; empty level entries inherit it* |
| `multiomics_cell_level_batch_key` | `null` | `null` | `resolve_batch_levels:cell_level` | Template: *one column: scGLUE use_batch + the Z_comp Harmony pass* |
| `multiomics_sample_level_batch_col` | `null` | `null` | `resolve_batch_levels:sample_level` | Template: *RMD reference group (x modality), sample-level Harmony, autotune, DE* |
| `multiomics_celltype_col` | `"cell_type"` | `"cell_type"` | `cell_types_multiomics:cell_type_column`<br>`propagate_cell_type:celltype_col`<br>`run_autotune:celltype_col`<br>`compute_sample_embedding:celltype_col`<br>`sample_distance:cell_type_column`<br>`run_trajectory_gam_differential_gene_analysis:celltype_col`<br>`proportion_test:celltype_col` | Joint cell typing writes its labels here. With `multiomics_cell_type_cluster: false` the labels must already be in this column. |
| `multiomics_modality_col` | `"modality"` | `"modality"` | `build_sample_adata:modality_col`<br>`cell_types_multiomics:modality_column`<br>`run_autotune:modality_col`<br>`compute_sample_embedding:modality_col`<br>`visualize_multimodal_embedding:modality_col` |  |

### Multi-omics: general multiomics settings

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_verbose` | `true` | `true` | 13 functions: `multiomics_preparation`, `harmonize_xglue`, `cell_types_multiomics`, `propagate_cell_type`, `run_autotune`, `compute_sample_embedding`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `raisinfit`, `visualize_multimodal_embedding`, `run_dimension_association_analysis` |  |
| `multiomics_use_gpu` | `true` | `true` | `harmonize_xglue:use_gpu`<br>`compute_sample_embedding:use_gpu` | GPU runs only if the global `use_gpu` is also `true`. |
| `multiomics_random_state` | `42` | `42` | `harmonize_xglue:random_state` | Seeds the `harmonize_xglue` pass only; the multi-omics sample embedding uses its built-in seed 42. |

### Multi-omics: GLUE preprocessing parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_ensembl_release` | `98` | `98` | `multiomics_preparation:ensembl_release` | Only used when `multiomics_integration`. |
| `multiomics_species` | `"homo_sapiens"` | `"homo_sapiens"` | `multiomics_preparation:species` | Only used when `multiomics_integration`. |
| `multiomics_use_highly_variable` | `true` | `true` | `multiomics_preparation:use_highly_variable` | Only used when `multiomics_integration`. |
| `multiomics_n_top_genes` | `2000` | `2000` | `multiomics_preparation:n_top_genes` | Only used when `multiomics_integration`. |
| `multiomics_n_top_peaks` | `50000` | `50000` | `multiomics_preparation:n_top_peaks` | Only used when `multiomics_integration`. |
| `multiomics_atac_min_cells_floor` | `10` | `10` | `multiomics_preparation:atac_min_cells_floor` | Only used when `multiomics_integration`. |
| `multiomics_n_pca_comps` | `50` | `50` | `multiomics_preparation:n_pca_comps` | Only used when `multiomics_integration`. |
| `multiomics_n_lsi_comps` | `50` | `50` | `multiomics_preparation:n_lsi_comps` | Only used when `multiomics_integration`. |
| `multiomics_lsi_n_iter` | `15` | `15` | *nothing* | **Dead key:** `multiomics_wrapper()` accepts `lsi_n_iter` but never reads it. |
| `multiomics_gtf_by` | `"gene_name"` | `"gene_name"` | `multiomics_preparation:gtf_by` | Only used when `multiomics_integration`. |
| `multiomics_flavor` | `"seurat_v3"` | `"seurat_v3"` | `multiomics_preparation:flavor` | Only used when `multiomics_integration`. |
| `multiomics_generate_umap` | `false` | `false` | `multiomics_preparation:generate_umap` | Only used when `multiomics_integration`. |
| `multiomics_compression` | `"gzip"` | `"gzip"` | *nothing* | **Dead key:** `multiomics_wrapper()` accepts `compression` but never reads it. |

### Multi-omics: GLUE training parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_consistency_threshold` | `0.05` | `0.05` | `multiomics_preparation:consistency_threshold` | Only used when `multiomics_integration`. |
| `multiomics_treat_sample_as_batch` | `false` | `false` | `multiomics_preparation:treat_sample_as_batch` | Only used when `multiomics_integration`. Template: *GLUE preserves per-sample variance (V2 default)* |
| `multiomics_save_prefix` | `"glue"` | `"glue"` | `multiomics_preparation:save_prefix` | Only used when `multiomics_integration`. |

### Multi-omics: scGLUE throughput knobs

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_glue_data_batch_size` | `1024` | **`128`** | `multiomics_preparation:data_batch_size` | Only used when `multiomics_integration`. |
| `multiomics_glue_max_epochs` | `null` | **`50`** | `multiomics_preparation:max_epochs` | Only used when `multiomics_integration`. |
| `multiomics_glue_dataloader_num_workers` | `4` | **`0`** | `multiomics_preparation:dataloader_num_workers` | Only used when `multiomics_integration`. |
| `multiomics_glue_dataloader_fetches_per_worker` | `8` | **`1`** | `multiomics_preparation:dataloader_fetches_per_worker` | Only used when `multiomics_integration`. |
| `multiomics_glue_array_shuffle_num_workers` | `4` | **`0`** | `multiomics_preparation:array_shuffle_num_workers` | Only used when `multiomics_integration`. |
| `multiomics_glue_graph_shuffle_num_workers` | `4` | **`0`** | `multiomics_preparation:graph_shuffle_num_workers` | Only used when `multiomics_integration`. |

### Multi-omics: cluster-vs-RMD split

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_harmonize_xglue_max_iter` | `50` | **`30`** | `harmonize_xglue:max_iter` |  |
| `multiomics_run_glue_twice_for_sample_removal` | `false` | `false` | `multiomics_preparation:run_second_glue_for_sample_removal` | Only used when `multiomics_integration`. |

### Multi-omics: cell-typing SNN label-transfer metric

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_metric` | `"cosine"` | `"cosine"` | `cell_types_multiomics:transfer_metric` | Only used when `multiomics_cell_type_cluster`. |

### Multi-omics: per-modality preprocess QC parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_rna_min_cells` | `500` | `500` | `multiomics_preparation:rna_min_cells` | Only used when `multiomics_integration`. |
| `multiomics_rna_min_genes` | `500` | `500` | `multiomics_preparation:rna_min_genes` | Only used when `multiomics_integration`. |
| `multiomics_rna_pct_mito_cutoff` | `20.0` | `20.0` | `multiomics_preparation:rna_pct_mito_cutoff` | Only used when `multiomics_integration`. |
| `multiomics_rna_exclude_genes` | `null` | `null` | `multiomics_preparation:rna_exclude_genes` | Only used when `multiomics_integration`. |
| `multiomics_atac_min_cells` | `1` | `1` | `multiomics_preparation:atac_min_cells` | Only used when `multiomics_integration`. |
| `multiomics_atac_min_features` | `2000` | `2000` | `multiomics_preparation:atac_min_features` | Only used when `multiomics_integration`. |
| `multiomics_atac_max_features` | `15000` | `15000` | `multiomics_preparation:atac_max_features` | Only used when `multiomics_integration`. |
| `multiomics_atac_min_cells_per_sample` | `1` | `1` | `multiomics_preparation:atac_min_cells_per_sample` | Only used when `multiomics_integration`. |
| `multiomics_atac_exclude_features` | `null` | `null` | `multiomics_preparation:atac_exclude_features` | Only used when `multiomics_integration`. |
| `multiomics_atac_doublet_detection` | `true` | `true` | `multiomics_preparation:atac_doublet_detection` | Only used when `multiomics_integration`. |
| `multiomics_atac_tfidf_scale_factor` | `10000.0` | `10000.0` | `multiomics_preparation:atac_tfidf_scale_factor` | Only used when `multiomics_integration`. |
| `multiomics_atac_log_transform` | `true` | `true` | `multiomics_preparation:atac_log_transform` | Only used when `multiomics_integration`. |

### Multi-omics: cell type clustering parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_existing_cell_types` | `false` | `false` | *nothing* | **Dead key:** `multiomics_wrapper()` accepts `existing_cell_types` but never reads it. |
| `multiomics_n_target_clusters` | `10` | **`8`** | *nothing* | **Dead key:** `multiomics_wrapper()` accepts `n_target_clusters` but never reads it. |
| `multiomics_cluster_resolution` | `0.8` | `0.8` | `cell_types_multiomics:cluster_resolution` | Only used when `multiomics_cell_type_cluster`. |
| `multiomics_use_rep_celltype` | `null` | `null` | `cell_types_multiomics:use_rep` | Only used when `multiomics_cell_type_cluster`. |
| `multiomics_markers` | `null` | `null` | *nothing* | **Dead key:** `multiomics_wrapper()` accepts `markers` but never reads it. |
| `multiomics_generate_umap_celltype` | `true` | `true` | `cell_types_multiomics:compute_umap` | Only used when `multiomics_cell_type_cluster`. |

### Multi-omics: visualization parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_plot_columns` | `null` | `null` | `multiomics_preparation:plot_columns` | Only used when `multiomics_integration`. |
| `multiomics_age_bin_size` | `null` | `null` | *nothing* | **Dead key:** only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call. |
| `multiomics_age_column` | `"age"` | `"age"` | *nothing* | **Dead key:** only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call. |
| `multiomics_plot_dendrogram_flag` | `true` | `true` | *nothing* | **Dead key:** only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call. |
| `multiomics_plot_cell_type_proportions_pca_flag` | `false` | `false` | *nothing* | **Dead key:** only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call. |
| `multiomics_plot_cell_type_expression_umap_flag` | `false` | `false` | *nothing* | **Dead key:** only used under `visualize_data`, which `wrapper()` hard-codes to `False` in this pipeline's `downstream_analysis()` call. |

### Multi-omics: sample embedding parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_derive_sample_embedding` | `true` | `true` | switch in multiomics_wrapper: run_autotune, compute_sample_embedding, save_embedding_to_h5ad, warn_if_stale_embedding |  |
| `multiomics_sample_embedding_medium_K` | `120` | `120` | `run_autotune:medium_K`<br>`compute_sample_embedding:medium_K` | Only used when `multiomics_derive_sample_embedding`. |
| `multiomics_sample_embedding_fine_K` | `300` | `300` | `run_autotune:fine_K`<br>`compute_sample_embedding:fine_K` | Only used when `multiomics_derive_sample_embedding`. |
| `multiomics_sample_embedding_rmd_dim` | `8` | `8` | `run_autotune:rmd_dim`<br>`compute_sample_embedding:rmd_dim_per_cluster` | Only used when `multiomics_derive_sample_embedding`. |
| `multiomics_sample_embedding_use_clr` | `false` | `false` | `compute_sample_embedding:use_clr` | Only used when `multiomics_derive_sample_embedding` and `not multiomics_autotune_enable`. Ignored when `multiomics_autotune_enable: true`. |
| `multiomics_sample_embedding_use_rmd` | `true` | `true` | `compute_sample_embedding:use_rmd` | Only used when `multiomics_derive_sample_embedding` and `not multiomics_autotune_enable`. Ignored when `multiomics_autotune_enable: true`. |
| `multiomics_sample_embedding_block_weights` | `null` | `null` | `compute_sample_embedding:block_weights` | Only used when `multiomics_derive_sample_embedding` and `not multiomics_autotune_enable`. Ignored when `multiomics_autotune_enable: true`. |
| `multiomics_sample_embedding_rmd_weight` | `"equal"` | `"equal"` | `compute_sample_embedding:rmd_weight` | Only used when `multiomics_derive_sample_embedding` and `not multiomics_autotune_enable`. Validated by the CLI: `"equal"` or a positive number. Ignored when `multiomics_autotune_enable: true`. Template: *RMD block weight α: "equal" = same energy as the three composition blocks together (α² = w_A1² + w_A2² + w_A3²); a number fixes α (0.60 = previous default)* |
| `multiomics_sample_embedding_pca_components` | `10` | `10` | `run_autotune:pca_components`<br>`compute_sample_embedding:pca_components` | Only used when `multiomics_derive_sample_embedding`. Capped at n_units - 1 (8 demo samples give 7 PCs). |
| `multiomics_sample_embedding_batch_method` | `"harmony"` | `"harmony"` | `run_autotune:batch_method`<br>`compute_sample_embedding:batch_method` | Only used when `multiomics_derive_sample_embedding`. `"harmony"`, `"linear"` or `"none"`; any other value is rejected by the CLI and by `compute_sample_embedding`. |

### Multi-omics: autotune parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_autotune_enable` | `false` | `false` | switch in multiomics_wrapper: run_autotune | Only used when `multiomics_derive_sample_embedding`. Replaces `compute_sample_embedding` with `run_autotune` (CPU only). Without a batch column and without `multiomics_autotune_grouping_col` no search runs and alpha stays "equal". |
| `multiomics_autotune_search` | `"bayesian"` | `"bayesian"` | `run_autotune:search` | Only used when `multiomics_autotune_enable` and `multiomics_derive_sample_embedding`. |
| `multiomics_autotune_scoring` | `"auto"` | `"auto"` | `run_autotune:scoring` | Only used when `multiomics_autotune_enable` and `multiomics_derive_sample_embedding`. |
| `multiomics_autotune_scope` | `"alpha_only"` | `"alpha_only"` | `run_autotune:scope` | Only used when `multiomics_autotune_enable` and `multiomics_derive_sample_embedding`. |
| `multiomics_autotune_alpha_bounds` | `[0.1, 100.0]` | `[0.1, 100.0]` | `run_autotune:alpha_bounds` | Only used when `multiomics_autotune_enable` and `multiomics_derive_sample_embedding`. |
| `multiomics_save_cell_adata_after_embedding` | `true` | `true` | `run_autotune:save_cell_adata`<br>`compute_sample_embedding:save_cell_adata`<br>switch in multiomics_wrapper: save_embedding_to_h5ad<br>switch in multiomics_wrapper: warn_if_stale_embedding | Only used when `multiomics_derive_sample_embedding`. Optional: may be left out of a config (default `true`). Template: *store .uns['X_DR_sample'] in the cell-level h5ad (updated in place; false = keep only the CSV)* |
| `multiomics_autotune_grouping_col` | `null` | `null` | `run_autotune:grouping_col` | Only used when `multiomics_autotune_enable` and `multiomics_derive_sample_embedding`. |
| `multiomics_autotune_tune_on_modality` | `null` | `null` | `run_autotune:tune_on_modality` | Only used when `multiomics_autotune_enable` and `multiomics_derive_sample_embedding`. |

### Multi-omics: trajectory analysis parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_trajectory_col` | `"sev.level"` | `"sev.level"` | `CCA_Call:trajectory_col`<br>`cca_pvalue_test:trajectory_col`<br>`run_trajectory_gam_differential_gene_analysis:anchor_col` | Only used when `multiomics_trajectory_analysis`. **Change for your data.** Default `sev.level`. Coerced to numeric for CCA: text labels become NaN and those samples are dropped. |
| `multiomics_trajectory_supervised` | `false` | **`true`** | `run_trajectory_gam_differential_gene_analysis:anchor_col`<br>switch in downstream_analysis: CCA_Call, cca_pvalue_test | Only used when `multiomics_trajectory_analysis`. |
| `multiomics_trajectory_visualization_label` | `null` | **`["sev.level"]`** | `TSCAN:grouping_columns` | Only used when `multiomics_trajectory_analysis` and `not multiomics_trajectory_supervised`. **Change for your data.** `null` becomes `["sev.level"]` (TSCAN plot colouring). |
| `multiomics_n_cca_pcs` | `2` | **`5`** | `CCA_Call:n_components`<br>`cca_pvalue_test:n_pcs` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_supervised`. Capped at the number of embedding PCs. Keep it well below the sample count: with PCs close to n_samples - 1 CCA fits any labelling perfectly. |
| `multiomics_cca_pvalue` | `false` | `false` | switch in downstream_analysis: cca_pvalue_test | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_supervised`. Permutation test with the same number of PCs as the observed score. When the PCs approach n_samples - 1 every permuted score is ~1: the step warns and the p-value is meaningless. |
| `multiomics_tscan_origin` | `null` | `null` | `TSCAN:origin` | Only used when `multiomics_trajectory_analysis` and `not multiomics_trajectory_supervised`. An integer: the index of the endpoint cluster to start from. |
| `multiomics_tscan_n_clusters` | `null` | `null` | `TSCAN:n_clusters` | Only used when `multiomics_trajectory_analysis` and `not multiomics_trajectory_supervised`. Template: *int; null lets TSCAN auto-select via BIC (2-20)* |
| `multiomics_tscan_pseudotime_mode` | `"rank"` | `"rank"` | `TSCAN:pseudotime_mode` | Only used when `multiomics_trajectory_analysis` and `not multiomics_trajectory_supervised`. `"rank"` or `"distance"`; anything else raises. Template: *"rank" or "distance"* |

### Multi-omics: sample distance parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_sample_distance_methods` | `null` | **`["cosine", "correlation"]`** | `sample_distance:method` | Only used when `multiomics_sample_distance_calculation`. SciPy metrics on the embedding (`euclidean`, `cosine`, `correlation`, … 12 in all) and `EMD`, `chi_square`, `jensen_shannon` on the cell-type composition; unknown names are skipped with a warning. |
| `multiomics_grouping_columns` | `null` | **`["sev.level"]`** | `sample_distance:grouping_columns` | Only used when `multiomics_sample_distance_calculation`. **Change for your data.** `null` becomes `["sev.level"]`. Every column must exist in the sample metadata, or the visualization step raises `ValueError: Missing required columns`. |
| `multiomics_summary_sample_csv_path` | `null` | `null` | `sample_distance:summary_csv_path` | Only used when `multiomics_sample_distance_calculation`. Written by the distance step and deleted after the trajectory step: always for the default path, and for a path you supply only if the run created it. |

### Multi-omics: trajectory DGE parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_fdr_threshold` | `0.05` | `0.05` | `run_trajectory_gam_differential_gene_analysis:fdr_threshold` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. |
| `multiomics_effect_size_threshold` | `1.0` | `1.0` | `run_trajectory_gam_differential_gene_analysis:effect_size_threshold` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. |
| `multiomics_top_n_genes` | `100` | `100` | `run_trajectory_gam_differential_gene_analysis:top_n_genes` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. |
| `multiomics_trajectory_diff_gene_covariate` | `null` | `null` | `run_trajectory_gam_differential_gene_analysis:covariate_columns` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. |
| `multiomics_trajectory_test` | `"gam"` | `"gam"` | `run_trajectory_gam_differential_gene_analysis:test` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. Template: *"gam" or "linear"; linear suits small cohorts* |
| `multiomics_num_splines` | `5` | `5` | `run_trajectory_gam_differential_gene_analysis:num_splines` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. |
| `multiomics_spline_order` | `3` | `3` | `run_trajectory_gam_differential_gene_analysis:spline_order` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. A gene is fitted only if it has a value in at least order + splines + 2 samples (after automatic reduction). With 8 samples use 2; order 3 fits nothing. |
| `multiomics_visualization_gene_list` | `null` | `null` | `run_trajectory_gam_differential_gene_analysis:visualization_gene_list` | Only used when `multiomics_trajectory_analysis` and `multiomics_trajectory_dge`. |

### Multi-omics: sample clustering parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_cluster_number` | `4` | **`3`** | `cluster:number_of_clusters` | Only used when `multiomics_sample_cluster`. |
| `multiomics_cluster_differential_gene_group_col` | `null` | **`"sev.level"`** | `proportion_test:group_col`<br>`raisinfit:group_col`<br>switch in downstream_analysis: proportion_test<br>switch in downstream_analysis: raisinfit, run_pairwise_tests | When set, the proportion test and RAISIN compare these groups **instead of** the k-means sample clusters (the demo sets `sev.level`). |

### Multi-omics: embedding visualization parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_color_col` | `null` | `null` | `visualize_multimodal_embedding:color_col`<br>switch in downstream_analysis: visualize_multimodal_embedding | Colour column of the modality-highlight plots (see `multiomics_target_modality`). |
| `multiomics_visualization_grouping_column` | `null` | `null` | `visualize_multimodal_embedding:visualization_grouping_column` | Without it the multi-omics embedding plots are uncoloured. |
| `multiomics_target_modality` | `"ATAC"` | `"ATAC"` | `visualize_multimodal_embedding:target_modality`<br>switch in downstream_analysis: visualize_multimodal_embedding | When set, the embedding plots also include modality-highlight plots `<modality>_<column>_{expression,proportion}.png`, coloured by `multiomics_color_col` or else the first grouping column. `null` skips them. |
| `multiomics_sample_embedding_key` | `"X_DR_sample"` | `"X_DR_sample"` | `visualize_multimodal_embedding:expression_key`<br>`visualize_multimodal_embedding:proportion_key` |  |
| `multiomics_figsize` | `[20, 8]` | `[20, 8]` | `visualize_multimodal_embedding:figsize` |  |
| `multiomics_point_size` | `60` | `60` | `visualize_multimodal_embedding:point_size` |  |
| `multiomics_alpha` | `0.8` | `0.8` | `visualize_multimodal_embedding:alpha` |  |
| `multiomics_colormap` | `"viridis"` | `"viridis"` | `visualize_multimodal_embedding:colormap` |  |
| `multiomics_show_sample_names` | `false` | `false` | `visualize_multimodal_embedding:show_sample_names` |  |
| `multiomics_force_data_type` | `null` | `null` | `visualize_multimodal_embedding:force_data_type` |  |

### Multi-omics: dimension association analysis parameters

| Key | Default | Demo | Consumed by | Notes |
|---|---|---|---|---|
| `multiomics_association_continuous_cols` | `null` | `null` | `run_dimension_association_analysis:continuous_cols` | Only used when `multiomics_dimension_association_analysis`. |
| `multiomics_association_categorical_cols` | `null` | `null` | `run_dimension_association_analysis:categorical_cols` | Only used when `multiomics_dimension_association_analysis`. |
| `multiomics_association_n_permutations` | `999` | **`99`** | `run_dimension_association_analysis:n_permutations` | Only used when `multiomics_dimension_association_analysis`. |

