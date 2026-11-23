# Downstream analysis

Every downstream function reads the sample embedding `uns['X_DR_sample']` (samples × PCs) that the [RNA tutorial](../rna.md) or the [Quick start](../quickstart.md) stores in the cell-level AnnData. This page runs each analysis once on the 8-sample COVID demo (4 healthy, `sev.level = 1`; 4 severe, `sev.level = 4`). The calls are the same for ATAC and multi-omics; only the input file changes.

!!! warning "8 samples is enough to run everything, not to conclude anything"
    With 8 samples, p-values sit at their floor, CCA fits almost any labels and no gene passes FDR. Use the demo to learn the calls and the files. The limits are listed under [Small cohorts](../downstream_details.md#small-cohorts).

## Setup

Run from the folder that holds `sampledisco_demo_output/`, in one Python session for the whole page:

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

OUT = "sampledisco_demo_output/rna"
adata = ad.read_h5ad(f"{OUT}/preprocess/adata_preprocessed.h5ad")  # cells x genes, with uns['X_DR_sample']
pseudo_adata = build_sample_adata(adata, sample_col="sample")      # one row per sample
print(adata.shape, pseudo_adata.shape)
print(pseudo_adata.obs["sev.level"].to_dict())
```

```text
[build_sample_adata] dropped 10 obs column(s) that vary within a sample (per-cell values): ['celltype', 'old', 'large', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'cell_type_original']
(29988, 12301) (8, 7)
{'CoV-3-Wilk': 4, 'CoV-34-Aruna': 4, 'CoV-36-Aruna': 4, 'CoV-4-Wilk': 4, 'HD-1-Wilk': 1, 'HD-20-Aruna': 1, 'HD-21-Aruna': 1, 'HD-4-Wilk': 1}
```

`adata` (cells) feeds the steps that need expression or cell types: EMD, trajectory DGE, the proportion test, RAISIN and the dendrogram. `pseudo_adata` (one row per sample, sample-level `obs` columns only) feeds everything else. The embedding has 7 PCs: SampleDisco caps the PC count at n_samples − 1. The whole page runs in 3 to 5 minutes and needs about 5 GB of RAM on 8 CPUs.

## What runs where

The config run (`sampledisco --config config.yaml`) runs these stages after the sample embedding. Flags carry the `rna_` prefix (`atac_` and `multiomics_` for the other branches), and folders are relative to `sampledisco_demo_output/rna/`. All stages go through [`downstream_analysis`](../../api/pipeline/downstream_analysis.md), and each stage's settings are in [Downstream options and caveats](../downstream_details.md).

| Analysis | Config flag | Python function | Output folder |
| --- | --- | --- | --- |
| [Sample distance](#sample-distance) | `rna_sample_distance_calculation` | `sample_distance` | `Sample_distance/<method>/` |
| [Dimension association](#dimension-association) | `rna_dimension_association_analysis` | `run_dimension_association_analysis` | `sample_association/` |
| [Trajectory: CCA](#trajectory) | `rna_trajectory_analysis`, `rna_trajectory_supervised: true`; `rna_cca_pvalue` | `CCA_Call`; `cca_pvalue_test` | `CCA/`; `CCA_test/` |
| [Trajectory: TSCAN](#trajectory) | `rna_trajectory_analysis`, `rna_trajectory_supervised: false` | `TSCAN` | `TSCAN/` |
| [Trajectory DGE](#trajectory-dge) | `rna_trajectory_dge` | `run_trajectory_gam_differential_gene_analysis` | `trajectoryDEG/` |
| [Sample clustering](#clustering) | `rna_sample_cluster` | `cluster` | `sample_cluster/` |
| [Proportion test](#clustering) | `rna_proportion_test` | `proportion_test` | `sample_cluster/proportion_test/` |
| [RAISIN DGE](#clustering) | `rna_cluster_dge` | `raisinfit`, `run_pairwise_tests` | `raisin_results/` |
| [Prediction, label transfer](#prediction) | none (Python only) | `predict_phenotype`, `transfer_labels` | none: they return tables |
| [Visualization](#visualization) | `rna_visualize_data` | `visualization` | `visualization/` |
| [Autotune α](#autotune) | `rna_autotune_enable` | `run_autotune` | `sample_embedding/` |

## Sample distance {#sample-distance}

`sample_distance` writes a sample × sample distance matrix, a clustered heatmap and a group-separation score (mean distance between groups / within groups; above 1, the groups separate). Cosine works on the embedding; EMD works on cell-type proportions, so it also needs the cells.

```python
from sampledisco.sample_distance.sample_distance import sample_distance

sample_distance(adata=pseudo_adata, output_dir=f"{OUT}/Sample_distance", method="cosine",
                data_type="RNA", grouping_columns=["sev.level"])
sample_distance(adata=pseudo_adata, output_dir=f"{OUT}/Sample_distance", method="EMD",
                cell_adata=adata,               # cells: cell-type proportions and centroids
                pseudobulk_adata=pseudo_adata,  # samples: group labels for the score
                sample_column="sample", cell_type_column="cell_type",
                grouping_columns=["sev.level"])
```

**Output** (under `Sample_distance/`): `cosine/distance_matrix_sample_DR.csv`, `cosine/sample_distance_sample_DR_heatmap.pdf` and the score in `cosine/distance_check_results_sample_DR_cosine.txt`; `EMD_distance/distance_matrix_EMD.csv` and its heatmap.

On the embedding, healthy samples sit closer to each other than to severe ones (cosine score 1.41). On cell-type proportions alone the groups barely separate (EMD score 1.02). All 15 methods are listed in [the details](../downstream_details.md#sample-distance).

![Cosine distances between the 8 samples](../../resource/downstream/ds_sample_distance_cosine.png)
<div class="figure-caption"><code>Sample_distance/cosine/sample_distance_sample_DR_heatmap.pdf</code>: cosine distances (scaled to a maximum of 1) with hierarchical clustering.</div>

## Dimension association {#dimension-association}

`run_dimension_association_analysis` reports how much of each PC's variance every sample-level variable explains (R², with a permutation p-value and an FDR).

```python
from sampledisco.sample_association.association import run_dimension_association_analysis

assoc = run_dimension_association_analysis(
    pseudo_adata=pseudo_adata, output_dir=f"{OUT}/sample_association",
    continuous_cols=[], categorical_cols=["sev.level", "batch"],  # the variables to test
    n_permutations=999, sample_col="sample", random_state=42,
)
table = assoc["results"]["X_DR_sample"]
print(table.sort_values("perm_p")[["variable", "component", "r2", "perm_p", "fdr"]].head(3))
```

```text
...
     variable component        r2  perm_p    fdr
13      batch       PC7  0.569886   0.022  0.308
2   sev.level       PC3  0.398277   0.119  0.574
0   sev.level       PC1  0.399612   0.123  0.574
```

**Output** (under `sample_association/`): `variance_explained_sample.csv` (one row per variable × PC), `figures/sample_variance_heatmap.pdf` and `figures/sample_top_associations.pdf`. On the demo, severity explains a large share of the leading PCs, but with 8 samples nothing survives the FDR.

## Trajectory: CCA, then TSCAN {#trajectory}

CCA (supervised) finds the direction in the first `K` PCs that correlates best with a numeric phenotype and gives every sample a pseudotime along it. `cca_pvalue_test` compares the score with label permutations; pass the same `K` as `n_pcs`.

```python
from sampledisco.sample_trajectory.CCA import CCA_Call
from sampledisco.sample_trajectory.CCA_test import cca_pvalue_test

K = 3                                     # PCs used by CCA; keep well below n_samples - 1
score, _, ptime, _ = CCA_Call(adata=pseudo_adata, output_dir=OUT,
                              trajectory_col="sev.level", n_components=K)
p = cca_pvalue_test(pseudo_adata=pseudo_adata, column="X_DR_sample", input_correlation=score,
                    output_directory=OUT, num_simulations=1000, trajectory_col="sev.level",
                    n_pcs=K)              # must equal n_components
print(f"CCA score on PC1-PC{K}: {score:.3f}, permutation p = {p:.3f}")
```

The output has the form:

```text
[CCA] Flipped pseudotime so it increases with 'sev.level' (Spearman rho before flip: ...)
P-value for observed correlation ...: ... [(1 + k) / (1000 + 1)]
CCA score on PC1-PC3: ..., permutation p = ...
```

The pseudotime `ptime` (`{sample: value}` in [0, 1]) is oriented to increase with `sev.level`. Expect a high CCA score and a small but not tiny p-value: with a 4 + 4 design an exact test cannot go below 1/35 ≈ 0.029, so the demo cannot show a strong effect. The exact numbers depend on the machine (see [run-to-run variation](../own_data.md#run-to-run-variation)).

TSCAN (unsupervised) orders the samples without a phenotype: it clusters them, joins the clusters with a minimum spanning tree and follows its longest path.

```python
from sampledisco.sample_trajectory.TSCAN import TSCAN

tscan = TSCAN(AnnData_sample=pseudo_adata, column="X_DR_sample", n_clusters=3, output_dir=OUT,
              grouping_columns=["sev.level"],  # colors a plot and orients the path
              origin=None, pseudotime_mode="rank")
print("cluster path:", tscan["main_path"], "origin:", tscan["origin"])  # [2, 0, 1], origin 2
```

**Output**: `CCA/pca_3d_cca_sample.pdf` (the samples on the PC pair best aligned with `sev.level`), `CCA/pseudotime_sample.csv`, `CCA_test/cca_pvalue_result_X_DR_sample.txt`, and `TSCAN/clusters_by_cluster_X_DR_sample.png`, `clusters_by_grouping_X_DR_sample.png` and `X_DR_sample_pseudotime.csv`. TSCAN's main path runs through clusters `[2, 0, 1]` (0-based), starting at the end with the lower mean `sev.level`; its pseudotime is in `tscan["pseudotime"]["main_path"]`.

![Samples on the chosen PC pair with the CCA direction](../../resource/downstream/ds_cca_sample.png)
<div class="figure-caption"><code>CCA/pca_3d_cca_sample.pdf</code>: the 8 samples on PC1 and PC3, colored by severity, with the CCA direction.</div>

## Trajectory DGE {#trajectory-dge}

`run_trajectory_gam_differential_gene_analysis` fits a GAM along a pseudotime to each cell type's pseudobulk expression and reports the genes that change.

```python
from sampledisco.sample_trajectory.trajectory_diff_gene import run_trajectory_gam_differential_gene_analysis

dge = run_trajectory_gam_differential_gene_analysis(
    adata,                          # cell level; the pseudobulk is built inside
    pseudotime_source=ptime,        # from CCA_Call; or tscan["pseudotime"]["main_path"]
    sample_col="sample", celltype_col="cell_type", batch_col="batch",
    num_splines=5, spline_order=2,  # the default spline_order=3 fits no gene on 8 samples
    fdr_threshold=0.05, effect_size_threshold=1.0, top_n_genes=100,
    anchor_col="sev.level",         # UP = rises with severity
    output_dir=f"{OUT}/trajectoryDEG", verbose=False,
)
print(f"{len(dge)} features fitted, {(dge['fdr'] < 0.05).sum()} with FDR < 0.05")  # 9882, 0
```

**Output** (under `trajectoryDEG/`, `<ts>` a timestamp): `gam_all_genes_<ts>.tsv` (one row per feature, named `<cell type> - <gene>`), `gam_pseudoDEGs_<ts>.tsv`, `gam_summary_<ts>.txt` and plots in `visualizations/`. The call warns that nothing is significant with 8 samples, as expected.

## Clustering, proportion test and RAISIN {#clustering}

`cluster` runs k-means on the embedding. The proportion test (CLR-transformed cell-type proportions, moderated t-test) and RAISIN (differential expression between groups of samples) then compare the clusters.

```python
import pandas as pd
from sampledisco.sample_clustering.cluster import cluster
from sampledisco.sample_clustering.proportion_test import proportion_test
from sampledisco.sample_clustering.RAISIN import raisinfit
from sampledisco.sample_clustering.RAISIN_TEST import run_pairwise_tests

sample_clusters, _ = cluster(pseudobulk_adata=pseudo_adata, output_dir=OUT,
                             number_of_clusters=2, random_state=0)  # keep every cluster >= 2 samples
print(pd.Series(sample_clusters).value_counts().to_dict())

props = proportion_test(adata=adata, sample_col="sample", sample_to_clade=sample_clusters,
                        celltype_col="cell_type",   # the default is "celltype"
                        output_dir=f"{OUT}/sample_cluster/proportion_test")

fit = raisinfit(adata=adata, sample_col="sample", sample_to_clade=sample_clusters,
                testtype="unpaired", n_jobs=8, verbose=False)
n_sig, tables = run_pairwise_tests(fit=fit, output_dir=f"{OUT}/raisin_results", verbose=False)
print(tables["0_vs_1"].sort_values("pvalue").head(3))
```

```text
...
{0: 6, 1: 2}
...
        Foldchange      stat    pvalue       FDR
RPS14     1.205664  4.077699  0.000045  0.559493
PNISR     0.608640  3.790880  0.000150  0.763767
S100A8   -1.885200 -3.736945  0.000186  0.763767
```

**Output**: `sample_cluster/kmeans_clusters_sample.csv` and `kmeans_sample_embedding.png`; `sample_cluster/proportion_test/proportion_test_0_vs_1.csv` (`logFC`, `p_value`, `FDR` per cell type) and proportion heatmaps; `raisin_results/0_vs_1/raisin_results.csv` (per gene), volcano plots and `summary_plots/`.

K-means splits off two severe samples, and with 8 samples neither test finds anything at FDR < 0.05. To compare phenotype groups instead of clusters, pass `group_col="sev.level"`, as the demo config does ([details](../downstream_details.md#group-col)). RAISIN is the heaviest step: 1 to 3 minutes and up to about 5 GB of RAM on 8 CPUs.

## Prediction and label transfer {#prediction}

`predict_phenotype` predicts a sample label from the embedding with cross-validation; `transfer_labels` labels new samples from their nearest labelled neighbours. Both return tables and write no files.

```python
from sampledisco.sample_prediction import predict_phenotype, transfer_labels

pred, metrics = predict_phenotype(pseudo_adata, label_col="sev.level",
                                  n_splits=4,       # at most the size of the smallest class
                                  random_state=0)
print({k: metrics[k] for k in ("n_samples", "n_pcs", "auc", "balanced_accuracy")})

query = pseudo_adata.obs_names.isin(["CoV-3-Wilk", "HD-4-Wilk"])  # pretend these are unlabelled
table, tmetrics = transfer_labels(pseudo_adata, label_col="sev.level", query=query, n_neighbors=3)
print(table[table["query"]])
```

```text
{'n_samples': 8, 'n_pcs': 7, 'auc': 0.6875, 'balanced_accuracy': 0.375}
            label  query  predicted  confidence
sample
CoV-3-Wilk      4   True          1    0.666667
HD-4-Wilk       1   True          1    1.000000
```

`pred` holds the out-of-fold prediction of every sample, and `table` the transferred labels with their confidence. With 4 samples per class, these scores are noise. For donors with several samples, see [`group_col` and `score_by`](../downstream_details.md#prediction).

## Visualization {#visualization}

`visualization` draws a cell-type dendrogram and the samples on PC1 and PC2, colored by each grouping column.

```python
from sampledisco.visualization.visualization_other import visualization

visualization(AnnData_cell=adata, pseudobulk_anndata=pseudo_adata,
              output_dir=OUT,            # the figures go to OUT/visualization/
              grouping_columns=["sev.level"], plot_dendrogram_flag=True,
              plot_cell_type_proportions_pca_flag=True, plot_cell_type_expression_umap_flag=True)
```

**Output** (under `visualization/`): `cell_type_dendrogram.png`, `expression_embedding_sev.level.png` and `proportion_embedding_sev.level.png` (the same PC1–PC2 scatter in two colormaps).

## Autotune {#autotune}

Optional. `run_autotune` chooses α, the weight of the RMD block in the embedding (default `"equal"`), by scoring how well each candidate separates `grouping_col` and mixes `batch_col`.

```python
from sampledisco.parameter_selection.autotune import run_autotune

tuned = run_autotune(
    adata.copy(),                        # the copy keeps adata's embedding unchanged
    output_dir="sampledisco_demo_output/rna_autotune",
    sample_col="sample", celltype_col="cell_type",
    batch_col="batch", grouping_col="sev.level", seed=42,
)
print(tuned["best_params"], round(tuned["best_score"], 3))  # {'rmd_weight': 100.0} 0.438
```

**Output**: `rna_autotune/sample_embedding/sample_embedding.csv` and `autotune_record.txt`. Here α lands on the upper bound of the search, and the score barely moves above α ≈ 40 (0.437 to 0.438): on 8 samples, read it as a mild preference for a heavier RMD block, not an exact value. To use a tuned embedding, pass `tuned["adata"]` to `build_sample_adata`. Search strategies, scorers and the config keys are in [the details](../downstream_details.md#autotune).
