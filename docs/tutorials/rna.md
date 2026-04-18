# RNA pipeline tutorial

This tutorial walks through the scRNA-seq branch from raw counts to the sample-level pseudobulk embedding. Every step shows the call you would make from a notebook, the files it writes, and the figure you should expect. Parameter values follow the canonical [`config_covid_rna.yaml`](https://github.com/).

The pipeline ends at the sample embedding. Everything after that (sample distance, trajectory, differential genes, clustering, ...) is a downstream task — see [Downstream analysis](downstream/index.md).

!!! note "Imports"
    The code below assumes `genodistance` is installed. Until the package is published, use the equivalent absolute imports from `code/` (e.g. `from preparation.rna_preprocess_gpu import preprocess_linux`).

## Inputs

- `RNA.h5ad` — cell-level raw counts; `.obs` must carry a sample column (default `"sample"`).
- `sample_meta.csv` — one row per sample keyed by `sample`, with phenotype columns such as `sev.level`, `age`, `batch`.

Output lands under `output_dir/rna/`.

## 1. Preprocessing

Read counts, merge metadata, QC-filter cells and genes, select HVGs, compute PCA, and run Harmony integration. Two AnnData objects come out: one for clustering (batch-corrected) and one for sample-level differential analysis (minimally processed, raw counts preserved).

```python
from genodistance.preparation import preprocess_linux

adata_cluster, adata_sample = preprocess_linux(
    h5ad_path="/data/test_RNA.h5ad",
    sample_meta_path="/data/sample_meta.csv",
    output_dir="/results/rna",
    sample_column="sample",
    cell_level_batch_key=None,
    min_cells=500,
    min_genes=500,
    pct_mito_cutoff=20,
    num_cell_hvgs=2000,
    cell_embedding_num_PCs=20,
    num_harmony_iterations=30,
    verbose=True,
)
```

**Writes** → `/results/rna/preprocess/adata_cell.h5ad`, `/results/rna/preprocess/adata_sample.h5ad`.

## 2. Cell-type clustering

Leiden clustering on the Harmony-integrated embedding, with optional UMAP. The resulting labels are written to `adata_cluster.obs["cell_type"]` and propagated to `adata_sample`.

```python
from genodistance.preparation import cell_types_linux

adata_cluster, adata_sample = cell_types_linux(
    anndata_cell=adata_cluster,
    anndata_sample=adata_sample,
    leiden_cluster_resolution=0.99,
    n_target_clusters=None,
    umap=True,
    save=True,
    output_dir="/results/rna",
    verbose=True,
)
```

**Writes** → updated h5ad files and a UMAP PNG under `preprocess/`.

![UMAP colored by Leiden cell types](../resource/rna/preprocess_cell_type.png)
<div class="figure-caption">Step 2 — UMAP of single cells colored by Leiden cell type assignments.</div>

## 3. Sample embedding

Aggregate single cells into a pseudobulk AnnData and compute two dimension reductions:

- `X_DR_expression` — HVG-based PCA of cell-type-pooled expression (Harmony-corrected if a batch column is supplied).
- `X_DR_proportion` — PCA of per-sample cell-type proportions.

```python
from genodistance.sample_embedding import calculate_sample_embedding

pseudo_dict, pseudo_adata = calculate_sample_embedding(
    adata=adata_sample,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,
    output_dir="/results/rna",
    sample_hvg_number=2000,
    n_expression_components=10,
    n_proportion_components=10,
    harmony_for_proportion=True,
    use_gpu=True,
    atac=False,
    save=True,
)
```

**Writes** → `/results/rna/pseudobulk/pseudobulk_sample.h5ad` (with `X_DR_expression` and `X_DR_proportion` in `.obsm`) plus CSV exports under `embeddings/`.

---

From here, everything else — sample distance, CCA / TSCAN trajectory, trajectory DGE, sample clustering, proportion test, RAISIN cluster DGE, visualization, and the optional CCA-guided resolution search — runs off this `pseudo_adata`. Continue to the [Downstream analysis tutorials](downstream/index.md).
