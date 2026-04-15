# ATAC pipeline tutorial

The scATAC-seq pipeline mirrors the RNA pipeline — same downstream stack, same dual-embedding output — but the preprocessing and clustering steps switch to TF-IDF normalization and LSI dimension reduction. This tutorial shows the calls step by step using the values from [`config_covid_rna.yaml`](https://github.com/) (ATAC block).

## Inputs

- `ATAC.h5ad` — cell × peak counts; `.obs` must carry a `sample` column.
- `sample_meta.csv` — per-sample metadata including any phenotype of interest (e.g. `sev.level`).

Output lands under `output_dir/atac/`.

## 1. Preprocessing

TF-IDF normalization, LSI projection, optional doublet removal, and Harmony integration. ATAC typically uses a high feature count (50,000) and 50 LSI components.

```python
from genodistance.preparation import preprocess_linux  # ATAC version

adata_cluster, adata_sample = preprocess_linux(
    h5ad_path="/data/test_ATAC.h5ad",
    sample_meta_path="/data/sample_meta.csv",
    output_dir="/results/atac",
    sample_column="sample",
    cell_embedding_num_PCs=50,
    num_cell_hvfs=50000,
    min_cells=1,
    min_features=2000,
    max_features=15000,
    tfidf_scale_factor=1e4,
    log_transform=True,
    drop_first_lsi=True,
    doublet_detection=True,
    num_harmony_iterations=30,
    verbose=True,
)
```

**Writes** → `/results/atac/preprocess/adata_cell.h5ad`, `adata_sample.h5ad`.

## 2. Cell-type clustering

`cell_types_linux` auto-detects `X_lsi_harmony` in `.obsm` and switches to cosine-metric neighborhood graphs for ATAC.

```python
from genodistance.preparation import cell_types_linux

adata_cluster, adata_sample = cell_types_linux(
    anndata_cell=adata_cluster,
    anndata_sample=adata_sample,
    leiden_cluster_resolution=0.8,
    n_target_clusters=None,
    umap=False,
    save=True,
    output_dir="/results/atac",
)
```

**Writes** → updated h5ad files.

A hierarchical view of the resulting cell types helps sanity-check the granularity:

![Cell-type dendrogram](../resource/atac/cell_type_dendrogram.png)
<div class="figure-caption">Step 2 — Hierarchical dendrogram across Leiden-derived ATAC cell types (produced later by the visualization step, shown here for context).</div>

## 3. Sample embedding

Set `atac=True` so that pseudobulk aggregation and dimension reduction use LSI-style processing instead of log-normalized PCA.

```python
from genodistance.sample_embedding import calculate_sample_embedding

pseudo_dict, pseudo_adata = calculate_sample_embedding(
    adata=adata_sample,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="/results/atac",
    sample_hvg_number=50000,
    n_expression_components=30,
    n_proportion_components=10,
    harmony_for_proportion=True,
    use_gpu=True,
    atac=True,
    save=True,
)
```

**Writes** → `/results/atac/pseudobulk/pseudobulk_sample.h5ad` with the two embeddings in `.obsm`.

## 4. Sample distance

Same call pattern as RNA; no ATAC-specific flag is needed because the embeddings already carry the right signal.

```python
from genodistance.sample_distance import sample_distance

for method in ["cosine", "correlation"]:
    sample_distance(
        adata=pseudo_adata,
        output_dir="/results/atac",
        method=method,
        data_type="ATAC",
        grouping_columns=["sev.level"],
    )
```

**Writes** → `/results/atac/Sample_distance/{cosine,correlation}/`.

![Sample distance heatmap (expression, cosine)](../resource/atac/sample_distance_expression_DR_heatmap_cosine.png)
![Sample distance heatmap (proportion, cosine)](../resource/atac/sample_distance_proportion_DR_heatmap_cosine.png)
<div class="figure-caption">Step 4 — Cosine distance heatmaps on the ATAC expression and proportion embeddings.</div>

## 5. Supervised trajectory (CCA)

For ATAC, `n_cca_pcs=2` is usually enough because LSI-based embeddings concentrate signal in the leading components.

```python
from genodistance.sample_trajectory import CCA_Call

cca_results = CCA_Call(
    adata=pseudo_adata,
    output_dir="/results/atac",
    trajectory_col="sev.level",
    n_components=2,
)
```

**Writes** → `/results/atac/CCA/pca_2d_cca_{expression,proportion}.pdf`.

![CCA on expression embedding](../resource/atac/cca_expression.png)
![CCA on proportion embedding](../resource/atac/cca_proportion.png)
<div class="figure-caption">Step 5 — 2D CCA projections for the ATAC expression and proportion embeddings.</div>

## 6. Sample clustering

```python
from genodistance import cluster

expr_clusters, prop_clusters = cluster(
    pseudobulk_adata=pseudo_adata,
    output_dir="/results/atac",
    number_of_clusters=4,
)
```

**Writes** → `/results/atac/sample_cluster/kmeans_clusters_{expression,proportion}.csv`.

![K-means clusters on expression embedding](../resource/atac/kmeans_expression_embedding.png)
![K-means clusters on proportion embedding](../resource/atac/kmeans_proportion_embedding.png)
<div class="figure-caption">Step 6 — K-means sample clusters.</div>

## 7. Proportion test

Compare cell-type compositions across the sample clusters from step 6 (or across any phenotype group).

```python
from genodistance.sample_clustering import proportion_test

proportion_test(
    adata=adata_sample,
    sample_col="sample",
    sample_to_clade=expr_clusters,
    celltype_col="cell_type",
    output_dir="/results/atac/sample_cluster/expression",
)
```

**Writes** → `/results/atac/sample_cluster/expression/proportion_test/` (heatmaps + significance matrix PNG + CSVs).

![Cell-type proportion heatmap grouped by cell type](../resource/atac/proportion_heatmap_group_by_celltype.png)
<div class="figure-caption">Step 7 — Sample cell-type proportions laid out by cell type, grouped by sample cluster. Asterisks indicate significant differential abundance.</div>

## Notes and what we skipped

- **Trajectory DGE for ATAC** can use the same `run_trajectory_gam_differential_gene_analysis` function, but it's usually more meaningful on gene-activity signal (RNA-like). See the [Downstream tutorial](downstream.md) for the full call.
- **Unsupervised TSCAN** works identically to RNA — just pass `pseudo_adata` from step 3.
- **Cluster DGE with RAISIN** is also shared; see Downstream.

Continue to the [Downstream tutorial](downstream.md) for those modules.
