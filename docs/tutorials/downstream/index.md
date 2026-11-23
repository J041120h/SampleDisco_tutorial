# Downstream analysis

Every downstream module reads one object: the sample embedding `uns['X_DR_sample']`, a samples × PCs DataFrame that `compute_sample_embedding` writes into the cell-level AnnData. The functions are the same for RNA, ATAC and multi-omics; only the input AnnData changes.

These pages run on the RNA branch of the 8-sample COVID-19 demo (4 healthy, `sev.level = 1`; 4 severe, `sev.level = 4`). Each page shows the Python call, the files it writes and what came out of our run.

!!! warning "8 samples is enough to run everything, not to conclude anything"
    CCA scores, permutation p-values, trajectory DGE and prediction scores are degenerate or very noisy with 8 samples. Use the demo to learn the calls and the outputs. For results that mean something, see the 32-donor [Sound Life RNA tutorial](../soundlife_rna.md).

## Config keys, functions and outputs

The config-driven run (`sampledisco --config config.yaml`, see the [Quick start](../quickstart.md)) runs these stages in the order of the table, after the sample embedding. Each stage has an on/off flag, shown here with the `rna_` prefix; the ATAC and multi-omics branches use the same keys with `atac_` and `multiomics_`. Folders are relative to the branch output, `sampledisco_demo_output/rna/` in the demo.

| Stage | Config flag → main keys | Function | Output folder | Page |
| --- | --- | --- | --- | --- |
| Autotune α (replaces the plain embedding) | `rna_autotune_enable` → `rna_autotune_{search,scoring,scope,alpha_bounds,grouping_col}` | `run_autotune` | `sample_embedding/` | [Autotune](autotune.md) |
| Sample distance | `rna_sample_distance_calculation` → `rna_sample_distance_methods`, `rna_grouping_columns`, `rna_summary_sample_csv_path` | `sample_distance` | `Sample_distance/<method>/` | [Sample distance](sample_distance.md) |
| Trajectory, supervised (CCA) | `rna_trajectory_analysis` with `rna_trajectory_supervised: true` → `rna_trajectory_col`, `rna_n_cca_pcs` | `CCA_Call` | `CCA/` | [CCA](trajectory_cca.md) |
| CCA permutation p-value | `rna_cca_pvalue` (CCA only) | `cca_pvalue_test` | `CCA_test/` | [CCA](trajectory_cca.md) |
| Trajectory, unsupervised (TSCAN) | `rna_trajectory_analysis` with `rna_trajectory_supervised: false` → `rna_tscan_{n_clusters,origin,pseudotime_mode}`, `rna_trajectory_visualization_label` | `TSCAN` | `TSCAN/` | [TSCAN](trajectory_tscan.md) |
| Trajectory DGE (runs only inside the trajectory stage) | `rna_trajectory_dge` → `rna_num_splines`, `rna_spline_order`, `rna_fdr_threshold`, `rna_effect_size_threshold`, `rna_top_n_genes`, `rna_trajectory_diff_gene_covariate`, `rna_visualization_gene_list` | `run_trajectory_gam_differential_gene_analysis` | `trajectoryDEG/` | [Trajectory DGE](trajectory_dge.md) |
| Sample clustering (k-means) | `rna_sample_cluster` → `rna_cluster_number` | `cluster` | `sample_cluster/` | [Clustering](cluster.md) |
| Proportion test | `rna_proportion_test` → `rna_cluster_differential_gene_group_col` | `proportion_test` | `sample_cluster/proportion_test/` | [Proportion test](proportion_test.md) |
| RAISIN differential expression | `rna_cluster_dge` → `rna_cluster_differential_gene_group_col` | `raisinfit` + `run_pairwise_tests` | `raisin_results/` | [RAISIN](raisin_dge.md) |
| Visualization | `rna_visualize_data` → `rna_plot_dendrogram_flag`, `rna_plot_cell_type_proportions_pca_flag`, `rna_plot_cell_type_expression_umap_flag`, `rna_grouping_columns`, `rna_age_bin_size`, `rna_age_column` | `visualization` | `visualization/` | [Visualization](visualization.md) |
| Dimension association | `rna_dimension_association_analysis` → `rna_association_{continuous_cols,categorical_cols,n_permutations}` | `run_dimension_association_analysis` | `sample_association/` | [Association](sample_association.md) |
| Prediction, label transfer | not in the config | `predict_phenotype`, `transfer_labels` | none (they return tables) | [Prediction](prediction.md) |

Points that the flags alone don't tell you:

- **`rna_cluster_differential_gene_group_col` replaces the clusters.** When it is set, the proportion test and RAISIN compare the groups of that column and ignore the k-means clusters. The demo config sets it to `sev.level`, so its `sample_cluster/proportion_test/` compares healthy with severe, not cluster with cluster.
- **Several keys default to the COVID demo.** `rna_grouping_columns` and `rna_trajectory_visualization_label` default to `["sev.level"]` and `rna_trajectory_col` to `"sev.level"`. On your own data, set them to columns that exist, or the visualization stage raises an error and stops the run.
- **Small-sample settings.** The demo config sets `rna_n_cca_pcs: 5` and `rna_spline_order: 2`; the [CCA](trajectory_cca.md) and [trajectory DGE](trajectory_dge.md) pages explain both. Its `rna_cluster_number: 3` left two single-sample clusters on our embedding; the [clustering page](cluster.md) uses 2.
- **Rerunning only the downstream stages from the config:** set `rna_preprocessing`, `rna_cell_type_cluster` and `rna_derive_sample_embedding` to `false`. The run then reads `rna/preprocess/adata_preprocessed.h5ad` (or `rna_adata_cell_path`) and uses the `uns['X_DR_sample']` already stored in it.

## Choosing a trajectory method

| | CCA (supervised) | TSCAN (unsupervised) |
| --- | --- | --- |
| Use when | You have a phenotype to follow: severity, age, time point | You have no phenotype, or want to see whether one emerges |
| Needs | A **numeric** `trajectory_col` (non-numeric values become NaN and those samples are dropped) | Nothing beyond the embedding; `grouping_columns` only colors the plots and orients the path |
| Key settings | `n_components`: the first k PCs. Keep k well below n_samples − 1; with k ≥ n − 1 the fit is perfect for any labels | `n_clusters` (`None` = choose by BIC), `origin` (an endpoint cluster, 0-based integer), `pseudotime_mode` (`"rank"` or `"distance"`) |
| Significance | `cca_pvalue_test` with `n_pcs` equal to `n_components` | none |
| Pseudotime | Every sample, projected on the best-aligned pair of PCs, scaled to [0, 1] and (since 0.6.0) oriented to increase with `trajectory_col` | Samples on the main path of the cluster tree; samples on side branches get a branch pseudotime only |

Both return a `{sample: pseudotime}` dict that [trajectory DGE](trajectory_dge.md) accepts.

## Before you start

The pages pick up where the [RNA tutorial](../rna.md) or the RNA config run stops: `sampledisco_demo_output/rna/preprocess/adata_preprocessed.h5ad`, a cell-level AnnData with cell types in `obs['cell_type']` and, usually, the sample embedding in `uns['X_DR_sample']`. Run every page from the folder that holds `sampledisco_demo_output/`, and start each Python session with this block:

```python
import anndata as ad
from sampledisco.sample_embedding import compute_sample_embedding
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

OUT = "sampledisco_demo_output/rna"
adata = ad.read_h5ad(f"{OUT}/preprocess/adata_preprocessed.h5ad")   # cells x genes

if "X_DR_sample" not in adata.uns:
    # Only when the RNA run did not store the embedding in the h5ad.
    # batch_col=["batch"] matches rna_sample_level_batch_col in the demo config.
    adata = compute_sample_embedding(
        adata, OUT, sample_col="sample", celltype_col="cell_type",
        batch_col=["batch"], use_gpu=False, save=True,
    )

pseudo_adata = build_sample_adata(adata, sample_col="sample")   # one row per sample
print(adata.shape, pseudo_adata.shape)
print(pseudo_adata.obs[["sev.level", "batch"]])
```

```text
(29988, 12301) (8, 7)
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

- `adata` is the cell-level object. The trajectory DGE, proportion test, RAISIN, EMD distance and the dendrogram read it.
- `pseudo_adata` has one row per sample, with the embedding in `uns['X_DR_sample']` and `obsm['X_DR_sample']`. Every other function reads it.
- The embedding has 7 PCs, not the requested 10: SampleDisco caps the PC count at n_samples − 1.
- Loading the h5ad took 2 to 5 s and 1.3 GB of RAM. When the embedding is missing, computing it took 1.5 to 10 s on 8 CPUs, and saving it back into the h5ad up to half a minute more on our network file system; the result matched the stored embedding exactly.

!!! warning "Only some columns of `pseudo_adata.obs` are sample-level"
    `build_sample_adata` copies **every** column of `adata.obs` and keeps the value of the first cell of each sample. For columns that are constant within a sample (`sample`, `batch`, `sev.level` in the demo) that is the sample's value. For cell-level columns (`n_genes`, `total_counts`, `pct_counts_mt`, `celltype`, `cell_type`, ...) it is one arbitrary cell's value, and it means nothing at the sample level. Name the real sample-level columns explicitly wherever a function takes metadata columns, as the [association page](sample_association.md) does.

## The whole downstream stack from Python

`downstream_analysis` (in `sampledisco/wrapper/wrapper.py`; see its [API page](../../api/pipeline/downstream_analysis.md)) is the function the config run calls after the embedding. It takes the same settings as the config keys above, without the `rna_` prefix, and runs every enabled stage on an embedding you already have. Some defaults differ from the demo config; this call sets the ones that matter on 8 samples:

```python
from sampledisco.wrapper.wrapper import downstream_analysis

res = downstream_analysis(
    pseudo_adata=pseudo_adata,
    output_dir="sampledisco_demo_output/rna_downstream",   # a separate folder
    modality="rna",
    status_flags={"rna": {}},
    adata_cell=adata,
    adata_sample=adata,             # cell-level object for the proportion test and RAISIN
    batch_col="batch",
    grouping_columns=["sev.level"],
    trajectory_supervised=True,     # CCA; False = TSCAN (the default)
    trajectory_col="sev.level",
    n_cca_pcs=3,
    spline_order=2,                 # the default 3 fits no gene on 8 samples
    cluster_number=2,               # the default 4 (and the demo's 3) leave 1-sample clusters
    proportion_test=True,
    dimension_association_analysis=True,
    association_continuous_cols=[],
    association_categorical_cols=["sev.level", "batch"],
)
print(res["status_flags"]["rna"])
```

```text
{'sample_distance_calculation': True, 'trajectory_analysis': True, 'trajectory_dge': True, 'sample_cluster': True, 'proportion_test': True, 'visualization': True, 'dimension_association_analysis': True}
```

It writes the folders of the table under `sampledisco_demo_output/rna_downstream/` and returns `{"pseudo_adata", "status_flags"}`. RAISIN (`cluster_DGE=True`) and the CCA p-value (`cca_pvalue=True`) are off by default. On 8 CPUs this call took 1 to 2.5 minutes depending on the node, most of it the trajectory DGE, and the session peaked at 3.7 GB including the DGE worker processes.

## Pages

| Page | Function | What it does |
| --- | --- | --- |
| [Autotune](autotune.md) | `run_autotune` | Chooses the RMD weight α of the embedding |
| [Sample distance](sample_distance.md) | `sample_distance` | Sample × sample distances: 12 vector metrics on the embedding, 3 on cell-type proportions |
| [Dimension association](sample_association.md) | `run_dimension_association_analysis` | Variance of each PC explained by each sample-level variable |
| [Trajectory: CCA](trajectory_cca.md) | `CCA_Call`, `cca_pvalue_test` | Pseudotime along a phenotype, with a permutation p-value |
| [Trajectory: TSCAN](trajectory_tscan.md) | `TSCAN` | Pseudotime from a cluster tree, without a phenotype |
| [Trajectory DGE](trajectory_dge.md) | `run_trajectory_gam_differential_gene_analysis` | Genes that change along a pseudotime (per cell type, GAM) |
| [Sample clustering](cluster.md) | `cluster` | K-means on the embedding |
| [Proportion test](proportion_test.md) | `proportion_test` | Cell-type proportions between groups (CLR + moderated t-test) |
| [RAISIN](raisin_dge.md) | `raisinfit`, `run_pairwise_tests` | Differential expression between groups of samples |
| [Prediction and label transfer](prediction.md) | `predict_phenotype`, `transfer_labels` | Cross-validated phenotype prediction; kNN label transfer |
| [Visualization](visualization.md) | `visualization` | Cell-type dendrogram and scatters of the embedding |

The pages write into the same folders as the config run (`sampledisco_demo_output/rna/...`) and overwrite them, so you can compare their outputs with the table above. The trajectory DGE page uses the pseudotime from the CCA page, and the proportion test and RAISIN pages use the clusters from the clustering page; every other page needs only the block above.
