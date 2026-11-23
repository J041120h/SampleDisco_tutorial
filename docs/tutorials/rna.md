# RNA

This page runs the scRNA-seq branch on the COVID demo, one function at a time: load the data, preprocess it, call cell types and build the sample embedding. The arguments are the demo-config values, so you get the RNA part of the [quick start](quickstart.md) run. It runs in about 3 minutes on 8 CPU cores and needs about 3 GB of RAM.

## 1. Load the data

Download the demo file and look at its samples; `preprocess` reads the file itself, so this step is only for a first look.

```bash
mkdir -p data
curl -fL -o data/test_RNA.h5ad "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
```

```python
import anndata as ad

adata = ad.read_h5ad("data/test_RNA.h5ad")
print(adata.shape)
print(adata.obs.groupby("sample", observed=True)[["batch", "sev.level"]].first())
```

**Output**

- `data/test_RNA.h5ad`: raw counts, 29,989 cells × 14,881 genes.
- 8 samples from two studies (`batch`: `Aruna`, `Wilk`) with a phenotype (`sev.level`: `1` = healthy, `4` = severe COVID-19). Each study has 2 samples of each level, so correcting `batch` keeps the phenotype.

## 2. Preprocess

`preprocess` filters cells and genes, normalizes, runs PCA and then Harmony twice: `Z_comp` removes `batch` and `sample` (used for cell typing), and `Z_rmd` removes only `batch`, keeping the sample differences the embedding measures.

```python
from sampledisco.preparation.rna_preprocess_cpu import preprocess

OUT = "sampledisco_demo_output/rna"

adata = preprocess(
    h5ad_path="data/test_RNA.h5ad",
    sample_meta_path=None,            # sample, batch and sev.level are already in .obs
    output_dir=OUT,
    sample_column="sample",
    cell_level_batch_key=["batch"],   # cell-level Harmony covariates
    sample_level_batch_key="batch",   # only checked to exist here; step 4 corrects it
    min_cells=100,
    min_genes=200,
    pct_mito_cutoff=20.0,
    exclude_genes=None,
    num_cell_hvgs=2000,
    cell_embedding_num_PCs=20,
    num_harmony_iterations=20,
    verbose=True,
)
print(adata.shape)
print(sorted(adata.obsm.keys()))
```

**Output**

- `adata`, also written to `sampledisco_demo_output/rna/preprocess/adata_preprocessed.h5ad`: 29,988 cells × 12,301 genes after QC, with `obsm['Z_comp']` and `obsm['Z_rmd']` (20 dimensions each).

## 3. Cell types

`cell_types` runs Leiden on `Z_comp`, merges the clusters down to 8 cell types and draws a UMAP.

```python
from sampledisco.preparation.cell_type_cpu import cell_types

adata = cell_types(
    anndata_cell=adata,
    cell_type_column="cell_type",
    existing_cell_types=False,        # cluster de novo; the file's 'celltype' labels are not used
    n_target_clusters=8,
    umap=True,
    save=True,
    output_dir=OUT,
    defined_output_path=f"{OUT}/preprocess/adata_preprocessed.h5ad",  # the default, written out
    leiden_cluster_resolution=0.8,
    cell_embedding_column=None,       # None = Z_comp
    cell_embedding_num_PCs=20,
    verbose=True,
)
print(adata.obs["cell_type"].value_counts().sort_index())
```

**Output**

- `adata.obs['cell_type']`: labels `"1"` to `"8"` (Leiden finds about 15–20 clusters, depending on the machine, and they are merged into 8), also saved in `adata_preprocessed.h5ad` with the UMAP.
- `sampledisco_demo_output/rna/preprocess/cell_type.csv` (`cell_id`, `cell_type`) and `umap_cell_type.png` (below).

![UMAP of the demo RNA cells colored by the 8 cell types](../resource/rna/umap_cell_type.png)
<div class="figure-caption">UMAP of the 29,988 cells, colored by the 8 cell types from this run.</div>

## 4. Sample embedding

`compute_sample_embedding` turns each sample's cells into one row: cell-type composition blocks plus an RMD block (how each cell type is shifted in that sample), reduced by PCA and corrected for `batch` at the sample level. Arguments not shown keep their defaults, which are also the demo-config values; [Sample embedding in depth](sample_embedding.md) explains the blocks and every option.

```python
from sampledisco.sample_embedding import compute_sample_embedding

adata = compute_sample_embedding(
    adata,
    OUT,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=["batch"],              # sample-level Harmony and RMD reference groups
    medium_K=120,
    fine_K=300,
    rmd_weight="equal",
    pca_components=10,
    seed=42,
)
emb = adata.uns["X_DR_sample"]
print(emb.shape)
```

**Output**

- `adata.uns['X_DR_sample']`: a DataFrame, 8 samples × 7 PCs (`pca_components` is capped at samples − 1), also saved in `adata_preprocessed.h5ad`.
- `sampledisco_demo_output/rna/sample_embedding/sample_embedding.csv` (the same table) and `sample_embedding_blocks.npz` (the blocks before weighting).

!!! note "Your numbers will differ a little"
    Reruns on the same kind of machine give identical results; on a different CPU type the cluster count and the embedding can move slightly. See [Run-to-run variation](own_data.md#run-to-run-variation).

Each argument above is an `rna_*` key of the demo config ([Configuration guide](configuration.md)); GPU versions of these calls are in [GPU acceleration](gpu.md).

**Next:** [Downstream analysis](downstream/index.md) runs distances, trajectories, clustering, association and prediction on `uns['X_DR_sample']`.
