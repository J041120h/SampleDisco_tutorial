# Downstream options and caveats

This page collects the options, limits and output details behind [Downstream analysis](downstream/index.md). Code blocks continue that page's Python session (its setup block, plus `sample_clusters` from its clustering block).

## Small cohorts {#small-cohorts}

Every analysis runs on 8 samples, but several hit hard limits there:

- **PCs ≤ n − 1.** The sample embedding keeps at most n_samples − 1 PCs (7 on the demo, whatever `pca_components` asks for), and every function caps its PC count at that.
- **CCA needs k well below n.** With k PCs, CCA is a regression of the phenotype on k predictors; at k ≥ n − 1 it fits any labelling, the score is 1.0 and every permuted score is about 1.0 too. On the demo, k = 2 to 7 gave scores of 0.67, 0.92, 0.95, 0.99, 1.00 and 1.00. Use k = 3 on 8 samples; the demo config's `rna_n_cca_pcs: 5` is valid but leaves fewer spare degrees of freedom.
- **Permutation p-values have a floor.** A 4 + 4 design has only 35 distinct splits of the labels, so an exact permutation p-value cannot go below 1/35 ≈ 0.029. A Monte-Carlo estimate below that (0.021 for the demo CCA) only reflects how many random permutations happened to draw the real split. The same floor applies to the association p-values, and an FDR over many tests cannot get far below it.
- **The GAM spline needs samples.** Trajectory DGE adapts the spline to the sample size, then skips every feature that is non-zero in fewer than spline order + number of splines + 2 samples:

    | `spline_order` | Order and splines used on 8 samples | Samples needed | On the demo |
    | --- | --- | --- | --- |
    | 3 (function default) | 3 and 4 | 9 | every feature is skipped; the table is empty |
    | 2 (demo config) | 2 and 3 | 7 | 9,882 of 15,688 features fitted |

    With 10 or more samples the default fits too.
- **Single-sample groups cannot be tested.** RAISIN estimates each group's between-sample variance from its samples; a group with one sample is listed in `fit["failgroup"]`, and `run_pairwise_tests` skips every pair that involves it (with a warning, a `FAILED` file in the pair's folder and a line in `raisin_summary.txt`). The proportion test skips pairs where a group has fewer than 2 samples, and lists them in `proportion_test_significant_summary.txt`. On the demo, k-means with k = 3 (the demo config) gave clusters of 6, 1 and 1 samples and both tests skipped all three pairs; k = 2 gave 6 and 2. Check the cluster sizes first.
- **Prediction folds.** `n_splits` must not exceed the smallest class (4 on the demo); the default of 5 raises `ValueError`.

## Sample distance: all 15 methods {#sample-distance}

`sample_distance` takes one `method` per call. Every matrix is divided by its largest entry, so values run from 0 to 1 and compare only within one method. An unknown name prints a warning and is skipped. `grouping_columns` names the column for the group-separation score (between-group / within-group mean distance).

| Methods | Input | Folder under `Sample_distance/` | Main files |
| --- | --- | --- | --- |
| `euclidean`, `sqeuclidean`, `minkowski`, `cityblock`, `chebyshev`, `cosine`, `correlation`, `canberra`, `braycurtis`, and the binary metrics `hamming`, `jaccard`, `matching` (of little use on PCs) | `pseudo_adata` (`scipy` `pdist` on `X_DR_sample`) | `<method>/` | `distance_matrix_sample_DR.csv`, `sample_distance_sample_DR_heatmap.pdf`, `distance_check_results_sample_DR_<method>.txt`, `distance_statistics_summary_<method>.csv` |
| `EMD`: earth mover's distance between cell-type proportions, with the ground distance between cell-type centroids in the first `n_pcs` (20) dimensions of `embedding_key` (`Z_comp`) | `cell_adata` (+ `pseudobulk_adata` for the score) | `EMD_distance/` | `distance_matrix_EMD.csv`, `cell_type_proportions.csv`, `cell_type_centroids.csv`, `distance_matrix_EMD_heatmap.pdf`, `distance_check_results_EMD_EMD.txt` |
| `chi_square`, `jensen_shannon`: distances between the proportion vectors | `cell_adata` (+ `pseudobulk_adata`) | `<method>/proportion_DR_distance/` | `distance_matrix_proportion_DR.csv`, `cell_type_abundances.pdf`, `sample_distance_proportion_DR_heatmap.pdf`, `distance_check_results_proportion_DR_<method>.txt` |

In the config, `rna_sample_distance_methods` lists the methods (the demo uses `["cosine", "correlation"]`). The config run passes the cell-level data, so the proportion methods work there too. With `rna_summary_sample_csv_path` set, each score is also added to that CSV.

## Dimension association

For each (variable, PC) pair, `run_dimension_association_analysis` fits a least-squares model: a continuous variable gives the squared Pearson correlation, a categorical one (one-hot) the ANOVA η². Each R² gets a permutation p-value (`n_permutations`) and a Benjamini–Hochberg FDR. If you leave `continuous_cols` and `categorical_cols` unset, every column of `pseudo_adata.obs` is classified and tested, so name the variables you care about. The config keys are `rna_association_continuous_cols`, `rna_association_categorical_cols` (both empty in the demo, which means auto) and `rna_association_n_permutations` (99 in the demo).

## Trajectories

### CCA

- **`trajectory_col` must be numeric.** It goes through `pd.to_numeric`; samples whose value cannot be converted are dropped from the score and the null.
- **Score, plot and pseudotime.** The score uses all of the first `n_components` PCs. The plot and the pseudotime use only the pair of PCs with the highest canonical correlation (`auto_select_best_2pc=True`). The pseudotime is scaled to [0, 1] and reversed when needed, so that it increases with `trajectory_col` (the log then shows a `Flipped pseudotime` line).
- **The p-value** refits the CCA on the first `n_pcs` PCs for `num_simulations` permutations and is (1 + k) / (N + 1), where k permuted scores reach the observed one. Set `n_pcs` equal to `n_components`. At the default `n_pcs=2` the null is built on fewer PCs than the score and the p-value sits at its minimum, 1/(N + 1), whatever the data.
- **Returns.** `CCA_Call` returns `(score, score, pseudotime, pseudotime)`; on failure it returns `(nan, nan, {}, {})` instead of raising, so check the score.
- **Config.** `rna_trajectory_supervised: true`, `rna_trajectory_col` (default `"sev.level"`), `rna_n_cca_pcs`; `rna_cca_pvalue: true` adds the test with `n_pcs = rna_n_cca_pcs`.

### TSCAN: `origin` and `pseudotime_mode`

TSCAN clusters the samples with a Gaussian mixture (`n_clusters=None` picks the number by BIC), joins the cluster centers with a minimum spanning tree and takes its longest path as the main trajectory.

- **`origin`** picks which end of the main path is time 0. It must be the 0-based index of one of the two endpoint clusters, or the call raises. With `origin=None`, TSCAN starts at the endpoint whose samples have the lower mean of the first numeric column in `grouping_columns`; without one, at the endpoint with the lower index. On the demo it started at the end with the lower mean `sev.level`.
- **`pseudotime_mode`**: `"rank"` gives the i-th of n ordered samples i/(n − 1); `"distance"` uses the cumulative distance along the path in the embedding, scaled to [0, 1].
- **Side branches.** Samples in clusters off the main path get `NaN` in `obs['tscan_pseudotime_main']` and a branch pseudotime in `tscan["pseudotime"]["branching_paths"]`.
- **Config.** `rna_trajectory_supervised: false`, with `rna_tscan_n_clusters`, `rna_tscan_origin` and `rna_tscan_pseudotime_mode`; `rna_trajectory_visualization_label` colors the plots.

Use CCA when you have a phenotype to follow (severity, age, time point), and TSCAN when you have none or want to see whether one emerges.

### Trajectory DGE

- **Pseudotime.** `pseudotime_source` takes a `{sample: pseudotime}` dict, a DataFrame, or the path of a CSV/TSV such as `CCA/pseudotime_sample.csv`.
- **Direction.** `UP` and `DOWN` follow the pseudotime. `anchor_col` names a numeric column and reverses the pseudotime when it runs against it; keep it for TSCAN, whose direction depends on `origin`. The config run sets it for TSCAN only.
- **Model.** Per cell type, the pseudobulk keeps the 2,000 most variable genes (`n_features_per_celltype`); with `batch_col` the batch is regressed out first. Each feature gets a GAM on the samples where it is non-zero, with a BH FDR over all features; the `top_n_genes` largest effects among the features with FDR < `fdr_threshold` are the pseudoDEGs. `n_jobs=-1` fits on every core.
- **Files.** Each run adds timestamped tables (`gam_all_genes_`, `gam_significant_`, `gam_pseudoDEGs_`, `gam_summary_`) instead of replacing them. The heatmap and cluster-pattern plots appear only with significant features; `visualization_gene_list` adds one plot per listed gene.
- **Config.** `rna_trajectory_dge` runs only inside the trajectory stage, with `rna_num_splines`, `rna_spline_order`, `rna_fdr_threshold`, `rna_effect_size_threshold`, `rna_top_n_genes`, `rna_trajectory_diff_gene_covariate` and `rna_visualization_gene_list`.

## Proportion test and RAISIN: `group_col` {#group-col}

Both tests take their groups from `sample_to_clade` (a `{sample: group}` dict such as the k-means clusters) or from `group_col` (a column of `adata.obs`). **When both are given, `group_col` wins.** The config run passes both: the clusters and `rna_cluster_differential_gene_group_col`. The demo config sets that key to `sev.level`, so its `proportion_test/` and `raisin_results/` compare healthy with severe samples, not clusters. This is that comparison from Python:

```python
by_severity = proportion_test(
    adata=adata, sample_col="sample", group_col="sev.level",  # overrides sample_to_clade
    celltype_col="cell_type", output_dir=f"{OUT}/sample_cluster/proportion_test_sev.level",
)
print(by_severity["1_vs_4"].head(3).round(4).to_string(index=False))

fit_sev = raisinfit(adata=adata, sample_col="sample", group_col="sev.level",
                    testtype="unpaired", n_jobs=8, verbose=False)
n_sig_sev, tables_sev = run_pairwise_tests(fit=fit_sev, output_dir=f"{OUT}/raisin_results_sev.level",
                                           verbose=False)
print(tables_sev["1_vs_4"].sort_values("pvalue").head(3))
```

```text
celltype   logFC  p_value    FDR
       2  1.6711   0.0083 0.0661
       4  1.1365   0.0552 0.2208
       6 -2.0099   0.0892 0.2380
...
      Foldchange      stat    pvalue       FDR
IGKC   -1.403446 -4.183352  0.000029  0.353337
HOPX    0.451075  3.459771  0.000541  1.000000
ENO1   -0.393033 -3.350979  0.000805  1.000000
```

Healthy samples have relatively more of cell type 2, and the immunoglobulin genes are higher in severe samples, but with 4 samples per group nothing passes FDR < 0.05.

Other points:

- **`celltype_col`** of `proportion_test` defaults to `"celltype"`, not `"cell_type"`. The demo h5ad also has a `celltype` column (the original authors' labels), so leaving it unset silently tests different cell types.
- **Proportion test output.** One `proportion_test_<g1>_vs_<g2>.csv` per pair (`celltype`, `logFC` as the difference in mean CLR, `p_value`, `FDR`), `proportion_test_significant_summary.txt`, proportion heatmaps (raw, CLR, z-score), a significance matrix, and box plots only for pairs with a cell type at FDR < 0.01.
- **RAISIN memory.** `raisinfit` makes a dense copy of the cells × genes matrix; on the demo the session peaked at 4 to 5 GB with `n_jobs=8`. With `batch_col`, ComBat corrects the expression first and needs more. The config run uses the first sample-level batch column and `n_jobs=-1` (every core the process may use, capped at `SLURM_CPUS_PER_TASK` in a SLURM job), and the demo config leaves `rna_cluster_dge` off.
- **RAISIN designs.** `testtype` also accepts `"paired"` (with `individual_col`), `"continuous"` and `"custom"` (with `custom_design`); `run_pairwise_tests(control_group=...)` tests every group against one control. Each pair gets `<g1>_vs_<g2>/raisin_results.csv` and volcano plots; `summary_plots/` holds the summary tables and, when some gene is significant, a pseudobulk heatmap and a dot plot.

## Prediction: `group_col` and `score_by` {#prediction}

`predict_phenotype` standardizes the first `n_pcs` PCs and fits an L2 logistic regression (`task="classification"`) or a ridge regression (`task="regression"`). The penalty is chosen by an inner 3-fold cross-validation inside each of the `n_splits` outer folds, and every sample is predicted once, by a model that did not see it.

- **`group_col`.** When a donor has several samples (repeat visits, or the RNA and ATAC units of one sample), pass the donor column. All samples of a group then stay in one fold, in the outer and the inner split; without it, a model can recognise the donor rather than the phenotype, and the scores are inflated.
- **`score_by="group"`** averages the out-of-fold probabilities (or predictions) and the labels within each `group_col` group before scoring, so each donor counts once. `pred` stays per sample either way.
- **Metrics.** `auc` (binary; the second class in sorted order is positive), `accuracy` and `balanced_accuracy`, or `r2` for regression, plus the penalty chosen in each fold. With few samples, the fold split alone moves an AUC by 0.1 or more; try several `random_state` values before comparing two embeddings.
- **`transfer_labels`** takes the majority label of the `n_neighbors` nearest labelled samples (Euclidean distance on the first `n_pcs` PCs; `weighted=True` weights votes by 1/distance). Every labelled sample is also labelled leave-one-out (`loo_accuracy`); with `group_col`, the other samples of its group are left out too. For multi-omics, pass `group_col="sample"`, since the RNA and ATAC units of one sample are each other's nearest neighbours.

Neither function is part of the config run.

## Visualization

- `visualization` adds a `visualization/` folder under `output_dir` itself, so pass the branch folder (`OUT`). Every grouping column must exist in `pseudo_adata.obs`, or the call raises.
- The two scatters (`plot_cell_type_proportions_pca_flag`, `plot_cell_type_expression_umap_flag`) show the same PC1–PC2 points in two colormaps; the names are historical and no PCA or UMAP is computed.
- The dendrogram clusters the cell-type means in `obsm['Z_comp']` (Ward linkage) and needs a cell-type column in the cell-level data.
- In the config, `rna_grouping_columns` (default `["sev.level"]`) colors the scatters; on your own data set it to a column that exists, or the visualization stage stops the run.

## Autotune {#autotune}

The embedding stacks three cell-type composition blocks and one RMD block. The composition weights are fixed by the block sizes; α (`rmd_weight`) weights the RMD block, and `"equal"` gives it the energy of the three composition blocks together. `run_autotune` searches α only (`scope="alpha_only"` is the only scope), rebuilding the embedding from blocks computed once.

| `scoring` | What is maximized |
| --- | --- |
| `"auto"` (default) = `"multi_metric_proxy"` | mean of the available proxies, each min-max scaled: grouping tracking (canonical correlation or R²) and SPS for `grouping_col`; iLISI and the negative batch silhouette for `batch_col` |
| `"ilisi_label"` | grouping tracking and batch iLISI only |
| `"cca"`, `"sps"`, `"cv_knn_severity"`, `"pseudotime_spearman"` | one grouping proxy |
| `"ilisi_batch"`, `"neg_asw_batch"` | one batch proxy |
| `"sev_minus_batch"` | grouping tracking − 0.5 × batch silhouette |

| `search` | Evaluations over `alpha_bounds` (log10 scale, default `(0.1, 100.0)`) |
| --- | --- |
| `"bayesian"` (default) | 5 log-spaced points including both bounds, then 10 by expected improvement |
| `"golden_section"` | 14 interior points, assuming one peak |
| `"grid"` | 15 log-spaced points including both bounds |

- **α on a bound.** `autotune_record.txt` then adds a `NOTE` line, and a range when several α score within 1% of the best. Widen `alpha_bounds` and search again; if the answer moves to the new bound while the score stays flat, there is no optimum: keep the default or any α on the plateau.
- **Circularity.** α is chosen to separate `grouping_col`; judge the tuned embedding on a label that took no part in the search.
- **Outputs.** `run_autotune` replaces `uns['X_DR_sample']` of the AnnData it gets, writes `<output_dir>/sample_embedding/sample_embedding.csv` and `autotune_record.txt` (`sample_embedding_tune-on-<modality>/` with `tune_on_modality`), and re-saves `<output_dir>/preprocess/adata_preprocessed.h5ad` when that file exists.

**From the config**, autotune runs instead of the plain embedding and needs:

- `rna_derive_sample_embedding: true` (with the stage off, autotune never runs) and `rna_autotune_enable: true`;
- something to score: `rna_autotune_grouping_col`, or `rna_sample_level_batch_col`. With neither, no search runs and α stays `"equal"`. With only the batch column (the demo config's default state), α is chosen for batch mixing alone.

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

The ATAC and multi-omics branches use the same keys with `atac_` and `multiomics_`; multi-omics adds `multiomics_autotune_tune_on_modality`. With autotune on, `*_sample_embedding_rmd_weight`, `*_sample_embedding_block_weights`, `*_sample_embedding_use_clr` and `*_sample_embedding_use_rmd` are ignored, and the embedding is built on the CPU even with `use_gpu: true`.

## Rerunning only the downstream stages

Set `rna_preprocessing`, `rna_cell_type_cluster` and `rna_derive_sample_embedding` to `false`. The run then reads `rna/preprocess/adata_preprocessed.h5ad` (or `rna_adata_cell_path`) and the `uns['X_DR_sample']` stored in it. From Python, [`downstream_analysis`](../api/pipeline/downstream_analysis.md) runs every enabled stage on an embedding you already have, with the same settings as the config keys minus the `rna_` prefix. `rna_trajectory_col` and `rna_grouping_columns` default to the demo's `sev.level`; set them to columns of your own data.
